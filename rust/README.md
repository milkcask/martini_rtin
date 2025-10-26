# martini_rtin

A Rust implementation of the RTIN (Right-Triangulated Irregular Network) algorithm for real-time terrain mesh generation.

[![Crates.io](https://img.shields.io/crates/v/martini_rtin.svg)](https://crates.io/crates/martini_rtin)
[![Documentation](https://docs.rs/martini_rtin/badge.svg)](https://docs.rs/martini_rtin)

Based on the paper ["Right-Triangulated Irregular Networks" by Will Evans et. al. (1997)](https://www.cs.ubc.ca/~will/papers/rtin.pdf) and inspired by [Mapbox's Martini library](https://github.com/mapbox/martini).

## Features

- Fast terrain mesh generation from height data
- Configurable level of detail based on error tolerance
- Memory-efficient hierarchical mesh representation
- No unsafe code

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
martini_rtin = "0.2.0"
```

## Example

```rust
use martini_rtin::Martini;

// Create a mesh generator for a 257x257 grid (2^8 + 1)
let martini = Martini::with_capacity(257);

// Generate terrain data (flat array of height values)
let terrain: Vec<f32> = (0..257*257).map(|i| {
    let x = i % 257;
    let y = i / 257;
    // Simple sine wave terrain
    ((x as f32 * 0.1).sin() + (y as f32 * 0.1).sin()) * 10.0
}).collect();

// Create a tile from the terrain data
let tile = martini.create_tile(terrain);

// Generate a mesh with maximum error of 1.0
let (vertices, triangles) = tile.get_mesh(1.0);

println!("Generated {} vertices and {} triangles", 
         vertices.len(), triangles.len() / 3);
```

## Algorithm

The RTIN algorithm works by:

1. Building a hierarchy of right triangles from the terrain grid
2. Computing approximation errors for each triangle level
3. Generating meshes by recursively subdividing triangles that exceed the error threshold

This approach allows for efficient level-of-detail mesh generation suitable for real-time applications.

## Grid Size Requirements

The grid size must be of the form 2^n + 1 (e.g., 3, 5, 9, 17, 33, 65, 129, 257, 513, 1025).

## Performance

The algorithm is designed for real-time use and can generate meshes from large terrain grids in milliseconds.

## License

ISC