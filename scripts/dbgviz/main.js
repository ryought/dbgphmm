/**
 * parameters
 */
var index = 0
var params = {
  target: '',
  time: 0,
}

var cy = window.cy = cytoscape({
  container: document.getElementById('cy'),
  style: [
    {
      selector: 'node',
      style: {
        shape: 'ellipse',
        label: function (e) { return e.data('label') },
        /*
        'width': function (e) { return e.data('sizes')[index] },
        'height': function (e) { return e.data('sizes')[index] },
        */
      }
    },
    {
      selector: 'edge',
      style: {
        // 'line-color': '#f92411',
        // 'width': '5px',
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


cy.on('click', 'node', function(e){
  const node = e.target;
  params.target = node.id();
})

var gui = new dat.GUI({name: 'My GUI'})
gui.add({name: 'sam'}, 'name')


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
const updateGraph = () => {
  /*
  cy.add({
    group: 'nodes',
    data: { id: i, label: i },
  })
  cy.add({
    group: 'edges',
    data: { id: i, source: i, target: 'n1', length: 1000 },
  })
  cy.add({
    group: 'edges',
    data: { id: i, source: i, target: 'n2', length: 1000 },
  })
  i += 1
  const layout = cy.layout({
    name: 'cola',
  })
  layout.run()
  */
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
gui.add(params, 'target').listen()


/**
 * layouts
 */
var layout
const start = () => {
  layout = cy.layout({
    name: 'cola',
    maxSimulationTime: 40000000,
  })
  layout.run()
}
const stop = () => {
  layout.stop()
}
gui.add({ start }, 'start')
gui.add({ stop }, 'stop')
