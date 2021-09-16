/**
 * parameters
 */
var index = 0
var params = {
  target: '',
  time: 0,
  depth: 10,
}
var gui = new dat.GUI({name: 'My GUI'})

var cy = window.cy = cytoscape({
  container: document.getElementById('cy'),
  style: [
    {
      selector: 'node',
      style: {
        shape: 'ellipse',
        // label: function (e) { return e.data('label') },
      }
    },
    {
      selector: 'edge',
      style: {
        label: function (e) { return e.data('label') },
        'curve-style': 'bezier',
        'target-arrow-shape': 'triangle',
      }
    }
  ],
  elements: fetch('/g.json')
    .then(function (res) { return res.json() })
    .then(function (json) { console.log('json', json); return json }),
  layout: {
    name: 'cola',
    maxSimulationTime: 1000,
  },
});


/**
 * node selection
 */
const select = (root, k) => {
  let depth = {}
  cy.elements().bfs({
    roots: cy.getElementById(root),
    visit: function(v, e, u, i, d) {
      depth[v.id()] = d
    },
    directed: false,
  })
  cy.nodes().style('display', 'none');
  cy.nodes().filter((v) => depth[v.id()] < k).neighborhood().style('display', '');
}
const unselect = () => {
  cy.nodes().style('display', '');
}
cy.on('click', 'node', function(e){
  const node = e.target;
  params.target = node.id();
  select(node.id(), params.depth);
})
gui.add(params, 'target').listen()
gui.add(params, 'depth', 1, 20, 1).listen()
gui.add({ unselect }, 'unselect')


/**
 * edge width history
 */
const updateWidth = () => {
  cy.edges()
    .style('width', function (e) {
      const widths = e.data('widths')
      if (params.time < widths.length) {
        return widths[params.time] * 5
      } else {
        return 0
      }
    })
}
const MAX_TIME = 10
gui.add(params, 'time', 0, MAX_TIME, 1)
  .listen()
  .onChange(() => {
    updateWidth()
  })
let timer = null;
const animate = () => {
  if (timer === null) {
    // start animation
    timer = setInterval(() => {
      updateWidth()
      params.time = (params.time + 1) % MAX_TIME
    }, 100)
  } else {
    // stop animation
    clearInterval(timer)
    timer = null
  }
}
gui.add({ animate }, 'animate')


/**
 * layouts
 */
let layout = null
const start = () => {
  if (layout === null) {
    const target = cy.elements().filter((element) => element.style('display') !== 'none')
    layout = target.layout({
      name: 'cola',
      maxSimulationTime: 40000000,
    })
    layout.run()
  }
}
const stop = () => {
  if (layout !== null) {
    layout.stop()
    layout = null
  }
}
gui.add({ start }, 'start')
gui.add({ stop }, 'stop')
