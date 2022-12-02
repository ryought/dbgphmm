var gui
var cy
const MAX_TIME = 10
var global_state = {
  // independent states
  selected_node: '',
  max_depth: 20,
  // synced states
  show_node_label: true,
  show_node_info: true,
  show_edge_label: true,
  edge_label_key: -1,
  edge_color_key: -1,
  edge_width_key: -1,
  node_label_key: -1,
  node_color_key: -1,
  node_size_key: -1,
  node_color_max: 1,
  use_history: true,
  time: 0,
  label: '',
}
var node_attrs = {}
var edge_attrs = {}

/**
 * parse attribute set from elements
 */
function parse_attrs(elements) {
  const node = elements
    .find((element) => element.group == 'nodes')
    .data.attrs.map((attr) => attr.type)
  const edge = elements
    .find((element) => element.group == 'edges')
    .data.attrs.map((attr) => attr.type)

  node_attrs = {
    'disabled': -1
  }
  for (let i = 0; i < node.length; i++) {
    node_attrs[`${i}: ${node[i]}`] = i
  }

  edge_attrs = {
    'disabled': -1
  }
  for (let i = 0; i < edge.length; i++) {
    edge_attrs[`${i}: ${edge[i]}`] = i
  }
}

/**
 * parse labels of history
 */
function parse_history_labels(elements) {
  const labels = elements
    .find((element) => element.group == 'history_labels')
  if (labels) {
    return labels.data.labels
  } else {
    return []
  }
}

/**
 *
 */
function get_attr(attrs, key, time) {
  if (key != null && key != -1) {
    const attr = attrs[key]
    if (attr.type === 'copy_nums') {
      return attr.value[time]
    } else {
      return attr.value
    }
  } else {
    return null
  }
}

/**
 * definition of color mapping x \in [x_min, x_max].
 */
function color(x, x_min, x_max) {
  const y = ((x - x_min) / (x_max - x_min)) * 255
  const z = Math.max(0, Math.min(y, 255))
  return `rgb(0, ${z}, ${z})`
}



/*
 * layouts
 */
// current ongoing layout
let layout = null;
function start_layout() {
  if (layout === null) {
    const target = cy.elements().filter((element) => element.style('display') !== 'none')
    layout = target.layout({
      name: 'cola',
      maxSimulationTime: 40000000,
      stop: () => {
        // remove layout when finished
        layout = null
      },
    })
    layout.run()
  }
}
function stop_layout() {
  if (layout !== null) {
    layout.stop()
    layout = null
  }
}



/*
 * node selection
 */
function select_node(root_node, max_depth) {
  let depth = {}
  const elem = cy.getElementById(root_node)
  if (elem) {
    cy.elements().bfs({
      roots: elem,
      visit: function(v, e, u, i, d) {
        depth[v.id()] = d
      },
      directed: false,
    })
    cy.nodes().style('display', 'none').style('text-background-opacity', 0.0);
    cy.nodes().filter((v) => depth[v.id()] < max_depth).neighborhood().style('display', '');
    elem.style('text-background-opacity', 0.5).style('text-background-color', 'green');
  }
}
function unselect_node() {
  cy.nodes().style('display', '');
}

/*
 * Json load
 */
function setup_file_input() {
  const input = document.getElementById('input');
  input.addEventListener('change', (event) => {
    const file = event.target.files[0];
    const reader = new FileReader();
    reader.addEventListener('loadend', (event) => {
      const text = event.target.result;
      const elements = JSON.parse(text);

      // run
      let labels = parse_history_labels(elements)
      parse_attrs(elements)
      init_cytoscape(elements)
      sync_states()
      init_controls(labels)
    });
    reader.readAsText(file);
  });
}
setup_file_input();
// function load_file_from_local() {
// 	return new Promise((resolve, reject) => {
//     const input = document.getElementById('input');
//     input.addEventListener('change', (event) => {
//       const file = event.target.files[0];
//       const reader = new FileReader();
//       reader.addEventListener('loadend', (event) => {
//         const text = event.target.result;
//         const json = JSON.parse(text);
//         resolve(json)
//       });
//       reader.readAsText(file);
//     });
//     input.click();
//   })
// }
// function load_file_from_remote() {
//   return fetch('/hoge.json')
//     .then((res) => res.json())
// }
function load_file() {
  gui.destroy()
  const input = document.getElementById('input');
  input.click()
}

/*
 * init functions
 */
function init_cytoscape(elements) {
  const container = document.getElementById('cy');
  container.classList.remove('hidden');
  cy = cytoscape({
    container: container,
    style: [
      {
        selector: 'node',
        style: {
          shape: 'ellipse',
          label: (e) => {
            if (e.scratch('show_label')) {
              const key = e.scratch('label_attr_key')
              const copy_num = e.data('copy_num')
              const info = e.scratch('show_info') ? e.data('info') : ''
              const label = `${e.data('label')} (x${copy_num})` || ''
              const use_history = e.scratch('use_history')
              const time = e.scratch('time')
              const history = use_history ? e.data('history')[time] : ''
              if (key != -1) {
                const attrs = e.data('attrs')
                const time = e.scratch('time')
                return `${label} (${get_attr(attrs, key, time)}) ${history} ${info}`
              } else {
                return `${label} ${history} ${info}`
              }
            } else {
              return ''
            }
          },
          'border-width': (e) => {
            const copy_num = e.data('copy_num')
            return Math.min(12, copy_num * 2)
          },
          'border-color': (e) => {
            // const copy_num = e.data('copy_num')
            // return color(copy_num + 1, 0, 5)
            return 'red'
          },
          'background-color': (e) => {
            const attrs = e.data('attrs')
            const use_history = e.scratch('use_history')
            const x_max = e.scratch('color_max')
            const time = e.scratch('time')
            if (use_history) {
              const history = e.data('history')
              return color(history[time], 0, x_max)
            } else {
              const key = e.scratch('color_attr_key')
              const x = get_attr(attrs, key, time)
              if (x != null)  {
                return color(x, 0, x_max)
              } else {
                return '#000'
              }
            }
          }
        }
      },
      {
        selector: 'edge',
        style: {
          label: (e) => {
            if (e.scratch('show_label')) {
              const key = e.scratch('label_attr_key')
              const label = e.data('label') || ''
              if (key != -1) {
                const attrs = e.data('attrs')
                const time = e.scratch('time')
                return `${label} (${get_attr(attrs, key, time)})`
              } else {
                return label
              }
            } else {
              return ''
            }
          },
          'line-color': (e) => {
            const attrs = e.data('attrs')
            const time = e.scratch('time')
            const key = e.scratch('color_attr_key')
            const x = get_attr(attrs, key, time)
            if (x != null) {
              return color(x, 0, 5)
            } else {
              return '#444'
            }
          },
          'width': (e) => {
            const attrs = e.data('attrs')
            const time = e.scratch('time')
            const key = e.scratch('width_attr_key')
            return get_attr(attrs, key, time) || 1
          },
          'curve-style': 'bezier',
          'target-arrow-shape': 'triangle',
        }
      }
    ],
    elements: elements.filter((e) => e.group == 'nodes' || e.group == 'edges'),
    layout: {
      name: 'random',
      maxSimulationTime: 1000,
    },
  })

  // add handlers
  cy.on('click', 'node', function (e) {
    const node = e.target
    global_state.selected_node = node.id()
    select_node(global_state.selected_node, global_state.max_depth)
  })
}

function sync_states() {
  cy.nodes().scratch('show_label', global_state.show_node_label)
  cy.nodes().scratch('show_info', global_state.show_node_info)
  cy.edges().scratch('show_label', global_state.show_edge_label)
  cy.edges().scratch('label_attr_key', global_state.edge_label_key)
  cy.edges().scratch('color_attr_key', global_state.edge_color_key)
  cy.edges().scratch('width_attr_key', global_state.edge_width_key)
  cy.nodes().scratch('label_attr_key', global_state.node_label_key)
  cy.nodes().scratch('color_attr_key', global_state.node_color_key)
  cy.nodes().scratch('size_attr_key', global_state.node_size_key)
  cy.nodes().scratch('color_max', global_state.node_color_max)
  cy.nodes().scratch('use_history', global_state.use_history)
  cy.elements().scratch('time', global_state.time)
}

function init_controls(history_labels) {
  gui = new dat.GUI({name: 'My GUI'})

  // [1] selection related
  const select = gui.addFolder('select')
  select.closed = false
  select.add({ load_file: load_file }, 'load_file')
  select.add(global_state, 'selected_node')
    .listen()
    .onChange(() => select_node(global_state.selected_node, global_state.max_depth))
  select.add(global_state, 'max_depth', 1, 40, 1)
  select.add({ unselect: unselect_node }, 'unselect')


  // [2] layout related
  const layout = gui.addFolder('layout')
  layout.closed = false
  layout.add({ start: start_layout }, 'start')
  layout.add({ stop: stop_layout }, 'stop')


  // [3] node/edge style related
  //
  const attr = gui.addFolder('attributes')
  attr.closed = false
  attr.add(global_state, 'show_node_label')
    .onChange((value) => cy.nodes().scratch('show_label', value))
  attr.add(global_state, 'show_node_info')
    .onChange((value) => cy.nodes().scratch('show_info', value))
  attr.add(global_state, 'show_edge_label')
    .onChange((value) => cy.edges().scratch('show_label', value))
  // edges
  attr.add(global_state, 'edge_label_key', edge_attrs)
    .onChange((value) => cy.edges().scratch('label_attr_key', value))
  attr.add(global_state, 'edge_color_key', edge_attrs)
    .onChange((value) => cy.edges().scratch('color_attr_key', value))
  attr.add(global_state, 'edge_width_key', edge_attrs)
    .onChange((value) => cy.edges().scratch('width_attr_key', value))
  // nodes
  attr.add(global_state, 'node_label_key', node_attrs)
    .onChange((value) => cy.nodes().scratch('label_attr_key', value))
  attr.add(global_state, 'node_color_key', node_attrs)
    .onChange((value) => cy.nodes().scratch('color_attr_key', value))
  attr.add(global_state, 'node_size_key', node_attrs)
    .onChange((value) => cy.nodes().scratch('size_attr_key', value))
  attr.add(global_state, 'node_color_max', 0, 30, 0.001)
    .onChange((value) => cy.nodes().scratch('color_max', value))


  // [4] animation related
  const animation = gui.addFolder('animation')
  animation.closed = false
  animation.add(global_state, 'use_history')
    .onChange((value) => cy.nodes().scratch('use_history', value))
  const n_history = history_labels.length
  const updateLabel = () => {
    if (history_labels.length > 0) {
      global_state.label = history_labels[global_state.time]
    } else {
      global_state.label = ''
    }
  }
  updateLabel()
  animation.add(global_state, 'time', 0, n_history - 1, 1)
    .onChange((time) => {
      cy.elements().scratch('time', time)
      updateLabel()
    })
  animation.add(global_state, 'label')
    .listen()
  // animation.add({ animate }, 'animate')
}




// gui.add(params, 'copy_num').listen()
// /**
//  * edge width history
//  */
// gui.add(params, 'time', 0, MAX_TIME, 1)
//   .listen()
//   .onChange(() => {
//     updateWidth()
//   })
// let timer = null;
// const animate = () => {
//   if (timer === null) {
//     // start animation
//     timer = setInterval(() => {
//       updateWidth()
//       params.time = (params.time + 1) % MAX_TIME
//     }, 100)
//   } else {
//     // stop animation
//     clearInterval(timer)
//     timer = null
//   }
// }


