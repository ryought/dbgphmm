var gui
var cy
const MAX_TIME = 10
var global_state = {
  // independent states
  selected_node: '',
  max_depth: 10,
  // synced states
  time: 0,
  show_node_label: true,
  show_edge_label: true,
  edge_label_key: null,
  edge_color_key: null,
  edge_width_key: null,
  node_label_key: null,
  node_color_key: null,
  node_size_key: null,
}
var node_attrs = {}
var edge_attrs = {}

function main() {
  fetch('/hoge.json')
    .then((res) => res.json())
    .then((elements) => {
      parse_attrs(elements)
      init_cytoscape(elements)
      sync_states()
      init_controls()
    })
}

main()

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
    'disabled': null
  }
  for (let i = 0; i < node.length; i++) {
    node_attrs[`${i}: ${node[i]}`] = i
  }

  edge_attrs = {
    'disabled': null
  }
  for (let i = 0; i < edge.length; i++) {
    edge_attrs[`${i}: ${edge[i]}`] = i
  }
}

function get_attr(attrs, key, time) {
  if (key != null) {
    const attr = attrs[key]
    if (attr.type === 'copy_nums') {
      return attr.value[t]
    } else {
      return attr.value
    }
  } else {
    return null
  }
}

function color(x, x_min, x_max) {
  const b = Math.floor(((x - x_min) / (x_max - x_min)) * 255)
  return `rgb(100, 100, ${b})`
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
    cy.nodes().style('display', 'none');
    cy.nodes().filter((v) => depth[v.id()] < max_depth).neighborhood().style('display', '');
  }
}
function unselect_node() {
  cy.nodes().style('display', '');
}



/*
 * init functions
 */
function init_cytoscape(elements) {
  cy = cytoscape({
    container: document.getElementById('cy'),
    style: [
      {
        selector: 'node',
        style: {
          shape: 'ellipse',
          label: (e) => {
            if (e.scratch('show_label')) {
              const key = e.scratch('label_attr_key')
              const label = e.data('label') || ''
              if (key != null) {
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
          'background-color': (e) => {
            const attrs = e.data('attrs')
            const time = e.scratch('time')
            const key = e.scratch('color_attr_key')
            const x = get_attr(attrs, key, time)
            if (x != null)  {
              return color(x, 0, 5)
            } else {
              return '#444'
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
              if (key != null) {
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
    elements: elements,
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
  cy.edges().scratch('show_label', global_state.show_edge_label)
  cy.edges().scratch('label_attr_key', global_state.edge_label_key)
  cy.edges().scratch('color_attr_key', global_state.edge_color_key)
  cy.edges().scratch('width_attr_key', global_state.edge_width_key)
  cy.nodes().scratch('label_attr_key', global_state.node_label_key)
  cy.nodes().scratch('color_attr_key', global_state.node_color_key)
  cy.nodes().scratch('size_attr_key', global_state.node_size_key)
  cy.elements().scratch('time', global_state.time)
}

function init_controls() {
  gui = new dat.GUI({name: 'My GUI'})

  // [1] selection related
  const select = gui.addFolder('select')
  select.closed = false
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


  // [4] animation related
  const animation = gui.addFolder('animation')
  animation.closed = false
  animation.add(global_state, 'time', 0, MAX_TIME, 1)
    .onChange((time) => {
      cy.elements().scratch('time', time)
    })
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


