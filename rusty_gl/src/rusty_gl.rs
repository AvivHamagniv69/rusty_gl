use std::f64::consts::PI;

use sdl2::rect::Point;
use sdl2::pixels::Color;

fn multiply_matricies(matrix_a: &[&[f64]]) {

}

#[derive(Copy, Clone)]
pub struct Point3d {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
impl Point3d {
    pub fn init(x: f64, y: f64, z: f64) -> Point3d {
        return Point3d { x: (x), y: (y), z: (z) };
    }
    
    pub fn project_point(point: Point3d, camera: &Camera) -> Point {
        let w = (2.0*camera.near_clipping_plane*camera.far_clipping_plane)/(camera.near_clipping_plane-camera.far_clipping_plane);
        let z = (point.z*(camera.far_clipping_plane/(camera.far_clipping_plane-camera.near_clipping_plane)) - ((camera.far_clipping_plane * camera.near_clipping_plane)/(camera.far_clipping_plane-camera.near_clipping_plane)));
        
        let x = ((camera.width as f64) / (camera.height as f64) * camera.fov * point.x)/z;
        let y = (camera.fov * point.y)/z;

        let new_x = x + (camera.width as f64/2.0);
        let new_y = y + (camera.height as f64/2.0);

        Point::new(new_x as i32, new_y as i32)
    }
    
    pub fn project_arr(points: &[Point3d], camera: &Camera) -> Vec<Point> {
        let mut vec: Vec<Point> = Vec::new();
        for point in points {
            vec.push(Point3d::project_point(Point3d { x: (point.x), y: (point.y), z: (point.z) }, camera))
        }
        return vec;
    }
    

}

#[derive(Debug)]
pub struct Buffer {
    pub points: Vec<Point>,
    pub indices: Vec<usize>,

    pub points_colors: Vec<Color>,
    pub lines_colors: Vec<Color>,
    pub triangles_colors: Vec<Color>,
}
impl Buffer {
    pub fn init() -> Buffer {
        let points: Vec<Point> = Vec::new();
        let indices: Vec<usize> = Vec::new();
        let points_color: Vec<Color> = Vec::new();
        let lines_color: Vec<Color> = Vec::new();
        let triangles_color: Vec<Color> = Vec::new();
        return Buffer { 
            points: (points),
            indices: (indices),
            points_colors: (points_color),
            lines_colors: (lines_color),
            triangles_colors: (triangles_color)
        }
    }

    pub fn load_indice_to_buffer(&mut self, indice: usize) {
        self.indices.push(indice);
    }

    pub fn load_indices_to_buffer(&mut self, indices: &[usize]) {
        for indice in indices {
            self.indices.push(*indice);
        }
    }

    pub fn load_vertex_to_buffer(&mut self, point: Point) {
        self.points.push(point);
    }
    
    pub fn load_vertecies_to_buffer(&mut self, points: &[Point]) {
        for point in points {
            self.points.push(*point);
        }
    }

    pub fn load_obj_to_buffer(&mut self, obj: &Obj3d, camera: &Camera) {
        for i in 0..obj.vertecies.len() {
            let projected_point = Point3d::project_point(obj.vertecies[i], camera);
            self.load_vertex_to_buffer(projected_point);
            self.points_colors.push(obj.color_points);
        }
        for i in 0..obj.indices.len() {
            self.load_indice_to_buffer(obj.indices[i]);
            self.lines_colors.push(obj.color_lines);
            self.triangles_colors.push(obj.color_triangles);
        }
    }
}

pub struct Obj3d {
    vertecies: Vec<Point3d>,
    indices: Vec<usize>,
    originPoint: Point3d,

    pub color_points: Color,
    pub color_lines: Color,
    pub color_triangles: Color,
}

impl Obj3d {
    pub fn init(points: &[Point3d], indices: &[usize], color_points: Color, color_lines: Color, color_triangles: Color) -> Obj3d {
        let points: Vec<Point3d> = points.to_vec();
        let indices: Vec<usize> = indices.to_vec();
        let originPoint = Obj3d::set_origin_point(&points);
        return Obj3d { vertecies: (points), indices: (indices), originPoint: (originPoint), color_points: (color_points), color_lines: (color_lines), color_triangles: (color_triangles) };
    }

    pub fn set_origin_point(points: &[Point3d]) -> Point3d {
        let mut origin_point: Point3d = Point3d::init(0.0, 0.0, 0.0);
        for point in points {
            origin_point.x += (*point).x;
            origin_point.y += (*point).y;
            origin_point.z += (*point).z;
        }
        return origin_point;
    }

    pub fn rotate_z_axis(&mut self, mut angle: f64) {
        angle *= PI/180.0;
        let sin_angle = f64::sin(angle);
        let cos_angle = f64::cos(angle);
        for i in 0..self.vertecies.len() {
            let y = self.vertecies[i].y - self.originPoint.y;
            let z = self.vertecies[i].z - self.originPoint.z;

            self.vertecies[i].y = (y * cos_angle - z * sin_angle) + self.originPoint.y;
            self.vertecies[i].z = (z * cos_angle + y * sin_angle) + self.originPoint.y;
        }
        self.originPoint = Obj3d::set_origin_point(&self.vertecies);
    }
}

pub struct Camera {
    pub origin_point: Point3d,
    pub width: u32,
    pub height: u32,

    fov: f64,
    fovRadians: f64,

    near_clipping_plane: f64,
    far_clipping_plane: f64,

    clip_matrix: [[f64; 4]; 4],
}

impl Camera {
    pub fn set_origin_point(&mut self, x: f64, y:f64, z: f64) {
        self.origin_point = Point3d {x: (x), y: (y), z: {z}};
    }

    pub fn set_fov(&mut self, fov: f64) {
        self.fov = fov;
        self.fovRadians = fov*PI/180.0;
    }

    pub fn set_near_clipping_plane(&mut self, near_clipping_plane: f64) {
        self.near_clipping_plane = 
        match near_clipping_plane {
            0.0..=1.0 => near_clipping_plane,
            _ => 0.1
        };
    }

    pub fn set_far_clipping_plane(&mut self, farClippingPlane: f64) {
        if (farClippingPlane >= 1000.0) {
            self.far_clipping_plane = farClippingPlane;
        } else {
            self.far_clipping_plane = 5000.0;
        }
    }

    pub fn init(origin_point: Point3d, width: u32, height: u32, mut fov: f64, nearClippingPlane: f64, farClippingPlane: f64) -> Camera {
        let temp_near_clipping_plane = 
        match nearClippingPlane {
            0.0..=1.0 => nearClippingPlane,
            _ => 0.1
        };

        let temp_far_clipping_plane = 
        if farClippingPlane >= 1000.0 {
            farClippingPlane
        } else {
            5000.0
        };

        fov = -(1.0 / f64::tan(fov/2.0));
        let clip_matrix: [[f64; 4]; 4] = 
        [
            [fov * (width as f64 /height as f64),0.0, 0.0, 0.0],
            [0.0, fov, 0.0, 0.0],
            [0.0, 0.0, (farClippingPlane + nearClippingPlane)/(farClippingPlane - nearClippingPlane), 1.0],
            [0.0, 0.0, (2.0*farClippingPlane * nearClippingPlane)/(farClippingPlane - nearClippingPlane), 0.0]
        ];

        return Camera {
            origin_point: (origin_point),
            width: (width),
            height: (height),
            fov: (fov),
            fovRadians: (fov*PI/180.0),
            near_clipping_plane: (temp_near_clipping_plane),
            far_clipping_plane: (temp_far_clipping_plane),
            clip_matrix: (clip_matrix),
        };
    }
}