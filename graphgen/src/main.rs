use rand::prelude::*;
use std::{collections::HashSet, fs::File, io::Write, cmp::{max, min}};

use clap::Parser;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[arg(short, long, default_value_t = 200000)]
    node: usize,

    /// Number of times to greet
    #[arg(short, long, default_value_t = 500000)]
    edge: usize,
}
#[derive(Debug, Default)]
struct Graph {
    node_cnt: usize,
    edge_list: Vec<(usize, usize, f64)>,
    exist: HashSet<(usize, usize)>,
}

impl Graph {
    fn write_to(&self, mut file: File) {
        file.write_fmt(format_args!(
            "{} {} {}\n",
            self.node_cnt,
            self.node_cnt,
            self.edge_list.len()
        ))
        .unwrap();
        for (f, t, w) in &self.edge_list {
            file.write_fmt(format_args!("{} {} {}\n", min(f, t) + 1, max(f, t) + 1, w))
                .unwrap();
        }
    }
    fn connected(&self, f: usize, t: usize) -> bool {
        self.exist.contains(&(f, t)) || self.exist.contains(&(t, f))
    }
    fn add_edge(&mut self, f: usize, t: usize, w: f64) -> bool {
        if !self.connected(f, t) && f != t {
            self.exist.insert((f, t));
            self.edge_list.push((f, t, w));
            true
        } else {
            false
        }
    }
}

fn main() {
    let Args { node, edge } = Args::parse();
    let mut g = Graph {node_cnt: node, edge_list: Vec::new(), exist: HashSet::new() };
    let mut rng = rand::thread_rng();
    for i in 1..node {
        let fa = rng.gen_range(0..i);
        g.add_edge(i, fa, rng.gen_range(1f64..5f64));
    }
    for _i in 0..edge - node + 1 {
        loop {
            let f = rng.gen_range(0..node);
            let t = rng.gen_range(0..node);
            if g.add_edge(f, t, rng.gen_range(1f64..5f64)) {
                break;
            }
        }
    }
    g.write_to(File::create(format!("g-{node}-{edge}.mtx")).unwrap());
}
