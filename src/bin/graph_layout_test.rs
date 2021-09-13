use forceatlas2::*;

const ITERATIONS: u32 = 2000;

fn main() {
    let nodes: usize = 4;
    let edges: Vec<(usize, usize)> = vec![(0, 1), (1, 2), (2, 3)];
    let mut layout = Layout::<f32>::from_graph(
        edges,
        Nodes::Degree(nodes),
        None,
        Settings {
            chunk_size: None,
            dimensions: 2,
            dissuade_hubs: false,
            ka: 0.01,
            kg: 0.001,
            kr: 0.002,
            lin_log: false,
            speed: 1.0,
            prevent_overlapping: None,
            strong_gravity: false,
        },
    );
    for _ in 0..ITERATIONS {
        layout.iteration();
    }
    for (i, pos) in layout.points.iter().enumerate() {
        println!("x={} y={}", pos[0], pos[1]);
    }
}
