var index = 0

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
    maxSimulationTime: 10000000,
  },
});

const layout = cy.layout({
  name: 'cola',
  randomize: true,
  // maxSimulationTime: 40000000,
  fit: true,
  infinity: true,
  // convergenceThreshold: 0.000000001,
  // edgeLength: function(edge){ return edge.data('length') / 10; },
})

cy.on('click', 'node', function(e){
  console.log('clicked', e)
})

var gui = new dat.GUI({name: 'My GUI'})
gui.add({name: 'sam'}, 'name')

layout.run()

var i = 0

const run = () => {
  /*
  cy.getElementById('n0').data('label', 'new label')
  cy.getElementById('e0').style('width', '20px')
  */
  index = (index + 1) % 3
  console.log('run')

  cy.nodes()
    .style('width', function (e) { return e.data('sizes')[index] })
    .style('height', function (e) { return e.data('sizes')[index] })

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
gui.add({ run }, 'run')

