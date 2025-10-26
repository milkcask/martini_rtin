//! Rust implementation of the RTIN (Right-Triangulated Irregular Network) algorithm
//! Based on the Martini library for terrain mesh generation

#![warn(missing_docs)]
#![doc = include_str!("../README.md")]

use glam::USizeVec2;

/// Martini - creates a mesh generator for a given grid size
pub struct Martini {
    grid_size: usize,
    num_triangles: usize,
    num_parent_triangles: usize,
    indices: Vec<usize>,
    coords: Vec<USizeVec2>,
}

impl Martini {
    /// Creates a new Martini mesh generator
    ///
    /// # Arguments
    /// * `grid_size` - Size of the terrain grid (must be 2^n+1, default 257)
    ///
    /// # Panics
    /// Panics if `grid_size` is not of the form 2^n+1
    pub fn with_capacity(grid_size: usize) -> Self {
        assert_valid_size(grid_size);

        let tile_size = grid_size - 1;
        let num_triangles = tile_size * tile_size * 2 - 2;
        let num_parent_triangles = num_triangles - tile_size * tile_size;

        let indices = vec![0usize; grid_size * grid_size];

        // coordinates for all possible triangles in an RTIN tile
        let mut coords = vec![USizeVec2::ZERO; num_triangles * 2];

        // Get triangle coordinates from its index in an implicit binary tree
        for i in 0..num_triangles {
            let mut id = i + 2;
            let (mut ax, mut ay, mut bx, mut by, mut cx, mut cy);

            if (id & 1) != 0 {
                // bottom-left triangle
                ax = 0;
                ay = 0;
                bx = tile_size;
                by = tile_size;
                cx = tile_size;
                cy = 0;
            } else {
                // top-right triangle
                ax = tile_size;
                ay = tile_size;
                bx = 0;
                by = 0;
                cx = 0;
                cy = tile_size;
            }

            while {
                id >>= 1;
                id > 1
            } {
                let mx = (ax + bx) >> 1;
                let my = (ay + by) >> 1;

                if (id & 1) != 0 {
                    // left half
                    bx = ax;
                    by = ay;
                    ax = cx;
                    ay = cy;
                } else {
                    // right half
                    ax = bx;
                    ay = by;
                    bx = cx;
                    by = cy;
                }
                cx = mx;
                cy = my;
            }

            let k = i * 2;
            coords[k] = USizeVec2::new(ax, ay);
            coords[k + 1] = USizeVec2::new(bx, by);
        }

        Self {
            grid_size,
            num_triangles,
            num_parent_triangles,
            indices,
            coords,
        }
    }

    /// Creates a tile from terrain data
    pub fn create_tile(&self, terrain: Vec<f32>) -> Tile<'_> {
        Tile::new(terrain, self)
    }
}

/// Tile - represents a terrain tile that can generate meshes at different error levels
pub struct Tile<'a> {
    terrain: Vec<f32>,
    martini: &'a Martini,
    errors: Vec<f32>,
}

impl<'a> Tile<'a> {
    /// Creates a new tile from terrain data
    ///
    /// # Arguments
    /// * `terrain` - Flat array of height values (size * size)
    /// * `martini` - Reference to the Martini mesh generator
    ///
    /// # Panics
    /// Panics if terrain data length doesn't match expected size
    pub fn new(terrain: Vec<f32>, martini: &'a Martini) -> Self {
        assert_terrain_len(&terrain, martini.grid_size);

        let mut tile = Self {
            errors: vec![0.0; terrain.len()],
            terrain,
            martini,
        };

        tile.update();
        tile
    }

    #[allow(clippy::many_single_char_names)]
    fn update(&mut self) {
        let num_triangles = self.martini.num_triangles;
        let num_parent_triangles = self.martini.num_parent_triangles;
        let coords = &self.martini.coords;
        let size = self.martini.grid_size;

        // Iterate over all possible triangles, starting from the smallest level
        for i in (0..num_triangles).rev() {
            let k = i * 2;
            let a = coords[k];
            let b = coords[k + 1];

            let m = USizeVec2::new((a.x + b.x) >> 1, (a.y + b.y) >> 1);
            let c = USizeVec2::new(m.x + m.y - a.y, m.y + a.x - m.x);

            // Calculate error in the middle of the long edge of the triangle
            let interpolated_height = f32::midpoint(
                self.terrain[a.y * size + a.x],
                self.terrain[b.y * size + b.x],
            );
            let middle_index = m.y * size + m.x;
            let middle_error = (interpolated_height - self.terrain[middle_index]).abs();

            self.errors[middle_index] = self.errors[middle_index].max(middle_error);

            if i < num_parent_triangles {
                // Bigger triangles; accumulate error with children
                let left_child_index = ((a.y + c.y) >> 1) * size + ((a.x + c.x) >> 1);
                let right_child_index = ((b.y + c.y) >> 1) * size + ((b.x + c.x) >> 1);
                self.errors[middle_index] = self.errors[middle_index]
                    .max(self.errors[left_child_index])
                    .max(self.errors[right_child_index]);
            }
        }
    }

    #[allow(clippy::too_many_lines)]
    /// Generates a mesh for the given maximum error
    ///
    /// # Arguments
    /// * `max_error` - Maximum allowed error for terrain approximation
    ///
    /// # Returns
    /// A tuple of (vertices, triangles) where:
    /// - vertices: flat array of `USizeVec2` vertex coordinates
    /// - triangles: flat array of triangle indices [a, b, c, a, b, c, ...]
    pub fn get_mesh(&self, max_error: f32) -> (Vec<USizeVec2>, Vec<usize>) {
        let size = self.martini.grid_size;
        let mut indices = self.martini.indices.clone();
        let mut num_vertices = 0usize;
        let mut num_triangles = 0usize;
        let max = size - 1;

        // Retrieve mesh in two stages that both traverse the error map:
        // Stage 1 - count elements: find used vertices and count triangles for allocation
        let count_elements = |indices: &mut Vec<usize>,
                              num_vertices: &mut usize,
                              num_triangles: &mut usize,
                              a: USizeVec2,
                              b: USizeVec2,
                              c: USizeVec2| {
            #[allow(clippy::too_many_arguments)]
            fn count_recursive(
                indices: &mut Vec<usize>,
                num_vertices: &mut usize,
                num_triangles: &mut usize,
                errors: &[f32],
                size: usize,
                max_error: f32,
                a: USizeVec2,
                b: USizeVec2,
                c: USizeVec2,
            ) {
                let m = (a + b) / 2;

                if (a.x.abs_diff(c.x) + a.y.abs_diff(c.y) > 1)
                    && errors[m.y * size + m.x] > max_error
                {
                    count_recursive(
                        indices,
                        num_vertices,
                        num_triangles,
                        errors,
                        size,
                        max_error,
                        c,
                        a,
                        m,
                    );
                    count_recursive(
                        indices,
                        num_vertices,
                        num_triangles,
                        errors,
                        size,
                        max_error,
                        b,
                        c,
                        m,
                    );
                } else {
                    let idx_a = a.y * size + a.x;
                    let idx_b = b.y * size + b.x;
                    let idx_c = c.y * size + c.x;

                    if indices[idx_a] == 0 {
                        *num_vertices += 1;
                        indices[idx_a] = *num_vertices;
                    }
                    if indices[idx_b] == 0 {
                        *num_vertices += 1;
                        indices[idx_b] = *num_vertices;
                    }
                    if indices[idx_c] == 0 {
                        *num_vertices += 1;
                        indices[idx_c] = *num_vertices;
                    }
                    *num_triangles += 1;
                }
            }

            count_recursive(
                indices,
                num_vertices,
                num_triangles,
                &self.errors,
                size,
                max_error,
                a,
                b,
                c,
            );
        };

        count_elements(
            &mut indices,
            &mut num_vertices,
            &mut num_triangles,
            USizeVec2::new(0, 0),
            USizeVec2::new(max, max),
            USizeVec2::new(max, 0),
        );
        count_elements(
            &mut indices,
            &mut num_vertices,
            &mut num_triangles,
            USizeVec2::new(max, max),
            USizeVec2::new(0, 0),
            USizeVec2::new(0, max),
        );

        let mut vertices = vec![USizeVec2::ZERO; num_vertices as usize];
        let mut triangles = vec![0usize; num_triangles * 3];
        let mut tri_index = 0usize;

        // Stage 2 - Process triangles: fill the vertex and triangle arrays
        let process_triangle = |vertices: &mut Vec<USizeVec2>,
                                triangles: &mut Vec<usize>,
                                tri_index: &mut usize,
                                a: USizeVec2,
                                b: USizeVec2,
                                c: USizeVec2| {
            #[allow(clippy::too_many_arguments)]
            fn process_recursive(
                vertices: &mut Vec<USizeVec2>,
                triangles: &mut Vec<usize>,
                tri_index: &mut usize,
                indices: &[usize],
                errors: &[f32],
                size: usize,
                max_error: f32,
                a: USizeVec2,
                b: USizeVec2,
                c: USizeVec2,
            ) {
                let m = (a + b) / 2;

                if (a.x.abs_diff(c.x) + a.y.abs_diff(c.y) > 1)
                    && errors[m.y * size + m.x] > max_error
                {
                    // Triangle doesn't approximate the surface well enough; drill down further
                    process_recursive(
                        vertices, triangles, tri_index, indices, errors, size, max_error, c, a, m,
                    );
                    process_recursive(
                        vertices, triangles, tri_index, indices, errors, size, max_error, b, c, m,
                    );
                } else {
                    // Add a triangle
                    let idx_a = indices[a.y * size + a.x] - 1;
                    let idx_b = indices[b.y * size + b.x] - 1;
                    let idx_c = indices[c.y * size + c.x] - 1;

                    vertices[idx_a] = a;
                    vertices[idx_b] = b;
                    vertices[idx_c] = c;

                    triangles[*tri_index] = idx_a;
                    triangles[*tri_index + 1] = idx_b;
                    triangles[*tri_index + 2] = idx_c;
                    *tri_index += 3;
                }
            }

            process_recursive(
                vertices,
                triangles,
                tri_index,
                &indices,
                &self.errors,
                size,
                max_error,
                a,
                b,
                c,
            );
        };

        process_triangle(
            &mut vertices,
            &mut triangles,
            &mut tri_index,
            USizeVec2::new(0, 0),
            USizeVec2::new(max, max),
            USizeVec2::new(max, 0),
        );
        process_triangle(
            &mut vertices,
            &mut triangles,
            &mut tri_index,
            USizeVec2::new(max, max),
            USizeVec2::new(0, 0),
            USizeVec2::new(0, max),
        );

        (vertices, triangles)
    }
}

/// Validates that the grid size is of the form 2^n+1
///
/// # Arguments
/// * `grid_size` - The grid size to validate
///
/// # Panics
/// Panics if the grid size is not of the form 2^n+1
///
/// # Examples
/// ```
/// use martini_rtin::assert_valid_size;
/// 
/// assert_valid_size(257); // Valid: 2^8 + 1
/// assert_valid_size(17);  // Valid: 2^4 + 1
/// ```
pub fn assert_valid_size(grid_size: usize) {
    assert!(
        grid_size > 2 && (grid_size - 1).is_power_of_two(),
        "Expected grid size to be 2^n+1, got {grid_size}",
    );
}

/// Validates that the terrain data length matches the expected grid size
///
/// # Arguments
/// * `terrain` - The terrain height data slice
/// * `grid_size` - The expected grid size (terrain should be grid_size × grid_size)
///
/// # Panics
/// Panics if the terrain length doesn't match grid_size²
///
/// # Examples
/// ```
/// use martini_rtin::assert_terrain_len;
/// 
/// let terrain = vec![0.0; 17 * 17];
/// assert_terrain_len(&terrain, 17); // Valid: 289 elements for 17×17 grid
/// ```
pub fn assert_terrain_len(terrain: &[f32], grid_size: usize) {
    assert_eq!(
        terrain.len(),
        grid_size * grid_size,
        "Expected terrain length to be {} ({} x {}) but got {}",
        grid_size * grid_size,
        grid_size,
        grid_size,
        terrain.len()
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_martini_creation() {
        let martini = Martini::with_capacity(257);
        assert_eq!(martini.grid_size, 257);
    }

    #[test]
    #[should_panic(expected = "Expected grid size to be 2^n+1")]
    fn test_invalid_grid_size() {
        Martini::with_capacity(256);
    }

    #[test]
    fn test_tile_creation() {
        let martini = Martini::with_capacity(17); // Small grid for testing
        let terrain = vec![0.0; 17 * 17];
        let _tile = martini.create_tile(terrain);
    }

    #[test]
    fn test_mesh_generation() {
        let martini = Martini::with_capacity(17);
        let terrain = vec![0.0; 17 * 17];
        let tile = martini.create_tile(terrain);
        let (vertices, triangles) = tile.get_mesh(0.0);

        // Should generate some mesh
        assert!(!vertices.is_empty());
        assert!(!triangles.is_empty());
        assert_eq!(triangles.len() % 3, 0); // Triangles come in groups of 3
    }
}
