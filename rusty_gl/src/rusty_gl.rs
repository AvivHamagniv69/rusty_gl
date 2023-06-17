use sdl2::render::Canvas;
use std::f64::consts::PI;
use std::str;
use std::fs::File;
use std::error::Error;
use std::io::BufReader;
use std::io::BufRead;
use std::time::Duration;

use sdl2;
use sdl2::pixels::Color;
use sdl2::rect::Point;
use sdl2::video::Window;

pub fn get_cross_product_f64(point1: &Point3, point2: &Point3) -> f64 {
    point1.x * point2.x + point1.y * point2.y + point1.z + point2.z
}

#[derive(Debug, Clone, Copy)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
impl Point3 {
    pub fn init(x: f64, y: f64, z: f64) -> Point3 {
        Point3 { x: (x), y: (y), z: (z)}
    }    

    fn project_point(point: Point3, camera: &Camera) -> Option<Point> {
        // a = forward
        // b = right
        // c = up
        let perspective = [
            1.0/((camera.width as f64/camera.height as f64)*f64::tan(camera.fov_radians as f64/2.0)),
            1.0/f64::tan(camera.fov_radians as f64/2.0),
            -((camera.far_clipping_plane+camera.near_clipping_plane)/(camera.far_clipping_plane-camera.near_clipping_plane)) as f64,
            -((2.0*camera.far_clipping_plane*camera.near_clipping_plane)/(camera.far_clipping_plane-camera.near_clipping_plane)) as f64
        ];

        let mut x: f64 = point.x-camera.origin_point.x;
        let mut y: f64 = point.y-camera.origin_point.y;
        let mut z: f64 = point.z-camera.origin_point.z;
        let mut w: f64 = -z;
        
        let temp_x = x;
        let temp_y = y;
        let temp_z = z;
        let temp_w = w;

        // checks if the point is behind the camera
        if z < 0.0 {
            return None;
        }

        x *= perspective[0];
        y *= perspective[1];
        z *= perspective[2] + perspective[3] * temp_z;

        x /= w;
        y /= w;

        x *= camera.fov_radians as f64;
        y *= camera.fov_radians as f64;

        x *= camera.width as f64;
        y *= camera.height as f64;

        x += camera.width as f64/2.0;
        y += camera.height as f64/2.0;

        return Some(Point::new(x as i32, y as i32));
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Triangle {
    pub point1: usize,
    pub point2: usize,
    pub point3: usize,

    origin_point: Point3,
    normal: Point3,

    pub color: Option<Color>,
}

pub struct Polygon {
    pub points: Vec<usize>,
    pub triangles: Vec<Triangle>,

    origin_point: Point3,
    normal: Point3,

    pub color: Option<Color>,
}

pub fn calc_origin_point_triangle(points: [Point3; 3]) -> Point3 {
    Point3 { 
        x: (points[0].x+points[1].x+points[2].x)/3.0,
        y: (points[0].y+points[1].y+points[2].y)/3.0,
        z: (points[0].z+points[1].z+points[2].z)/3.0
    }
}

pub fn normalize_vector(point: Point3) -> Point3 {
    let m = f64::sqrt(f64::powi(point.x, 2) + f64::powi(point.y, 2)+ f64::powi(point.z, 2));

    Point3 { 
        x: (point.x/m),
        y: (point.y/m),
        z: (point.z/m)
    }
}

pub fn calc_normal(points: [Point3; 3]) -> Point3 {
    let vu = Point3::init(points[1].x-points[0].x, points[1].y-points[0].y, points[1].z-points[0].z);
    let vv = Point3::init(points[2].x-points[0].x, points[2].y-points[0].y, points[2].z-points[0].z);

    let a = Point3 {
        x: (vu.y*vv.z) - (vu.z*vv.y),
        y: (vu.z*vv.x) - (vu.x*vv.z),
        z: (vu.x*vv.y) - (vu.y*vv.x)
    };
    normalize_vector(a)
}

pub trait MeshAndBuffer {
    fn sort_triangles_origin_point_z_distance_from_camera(&mut self, camera: &Camera, start: isize, end: isize) {}
    fn sort_triangles_in_mesh_origin_point_z_distance_from_camera(&mut self, camera: &Camera) {}

    fn sort_points_z_distance_from_camera(&mut self, camera: &Camera, start: isize, end: isize) {}
    fn sort_points_in_mesh_z_distance_from_camera(&mut self, camera: &Camera) {}
}

#[derive(Debug, Clone)]
pub struct Mesh {
    pub name: String,

    pub points: Vec<Point3>,
    pub triangles: Vec<Triangle>,

    pub origin_point: Point3,
    color: Option<Color>
}

impl MeshAndBuffer for Mesh {
    fn sort_triangles_origin_point_z_distance_from_camera(&mut self, camera: &Camera, start: isize, end: isize) {
        if end <= start {return;}

        let mut i = start as usize;
        let mut j: isize = i as isize-1;
        let pivot = end as usize;

        while i < pivot {
            let length_point_i_to_camera = camera.origin_point.z - self.triangles[i].origin_point.z;
            let length_point_pivot_to_camera = camera.origin_point.z - self.triangles[pivot].origin_point.z;

            if length_point_i_to_camera > length_point_pivot_to_camera {
                j += 1;
                self.triangles.swap(i, j as usize);
            }
            i += 1;
        }
        j += 1;
        self.triangles.swap(j as usize, pivot);

        self.sort_triangles_origin_point_z_distance_from_camera(camera, start, j -1);
        self.sort_triangles_origin_point_z_distance_from_camera(camera, j + 1, end);
    }
    fn sort_triangles_in_mesh_origin_point_z_distance_from_camera(&mut self, camera: &Camera) {
        if self.triangles.len() == 0 {return;}
        let a = self.triangles.len()-1;
        self.sort_triangles_origin_point_z_distance_from_camera(camera, 0, a as isize);
    }

    fn sort_points_z_distance_from_camera(&mut self, camera: &Camera, start: isize, end: isize) {
        if end <= start {return;}

        let mut i = start as usize;
        let mut j: isize = i as isize-1;
        let pivot = end as usize;

        while i < pivot {
            let length_point_i_to_camera = camera.origin_point.z - self.points[i].z;
            let length_point_pivot_to_camera = camera.origin_point.z - self.points[pivot].z;

            if length_point_i_to_camera > length_point_pivot_to_camera {
                j += 1;
                self.points.swap(i, j as usize);
            }
            i += 1;
        }
        j += 1;
        self.points.swap(j as usize, pivot);

        self.sort_points_z_distance_from_camera(camera, start, j -1);
        self.sort_points_z_distance_from_camera(camera, j + 1, end);
    }
    fn sort_points_in_mesh_z_distance_from_camera(&mut self, camera: &Camera) {
        let a = self.points.len()-1;
        self.sort_points_z_distance_from_camera(camera, 0, a as isize);
    }
}
impl Mesh {
    pub fn init(name: String, points: &[Point3], triangles: &[Triangle], color: Color) -> Mesh {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        for i in points {
            x += i.x;
            y += i.y;
            z += i.z;
        }
        let p = Point3 {x: (x/points.len() as f64), y: (y/points.len() as f64), z: (z/points.len() as f64)};
        Mesh {name: (name), points: (points.to_vec()), triangles: (triangles.to_vec()), origin_point: (p), color: (Some(color))}
    }

    pub fn set_origin_point(&mut self) {
        let mut p = Point3::init(0.0, 0.0, 0.0);
        for i in 0..self.points.len() {
            p.x += self.points[i].x;
            p.y += self.points[i].y;
            p.z += self.points[i].z;
        }
        p.x /= self.points.len() as f64;
        p.y /= self.points.len() as f64;
        p.z /= self.points.len() as f64;
        self.origin_point = p;
    }

    pub fn rotate_x(&mut self, mut angle: f64) {
        angle *= PI/180.0;
        self.set_origin_point();
        for p in 0..self.points.len() {
            self.points[p].y -= self.origin_point.y;
            self.points[p].z -= self.origin_point.z;

            let y = self.points[p].y;
            let z = self.points[p].z;

            self.points[p].y = y * f64::cos(angle) - z * f64::sin(angle);
            self.points[p].z = y * f64:: sin(angle) + z * f64::cos(angle);

            self.points[p].y += self.origin_point.y;
            self.points[p].z += self.origin_point.z;
        }
        self.set_origin_point();
        self.set_normals();
    }

    pub fn rotate_y(&mut self, mut angle: f64) {
        angle *= PI/180.0;
        self.set_origin_point();
        for p in 0..self.points.len() {
            self.points[p].x -= self.origin_point.x;
            self.points[p].z -= self.origin_point.z;

            let x = self.points[p].x;
            let z = self.points[p].z;

            self.points[p].x = x * f64::cos(angle) + z * f64::sin(angle);
            self.points[p].z = -(x * f64:: sin(angle)) + z * f64::cos(angle);

            self.points[p].x += self.origin_point.x;
            self.points[p].z += self.origin_point.z;
        }
        self.set_origin_point();
        self.set_normals();
    }

    pub fn rotate_z(&mut self, mut angle: f64) {
        angle *= PI/180.0;

        self.set_origin_point();
        for p in 0..self.points.len() {
            self.points[p].x -= self.origin_point.x;
            self.points[p].y -= self.origin_point.y;

            let x = self.points[p].x;
            let y = self.points[p].y;

            self.points[p].x = x * f64::cos(angle) - y * f64::sin(angle);
            self.points[p].y = x * f64:: sin(angle) + y * f64::cos(angle);

            self.points[p].x += self.origin_point.x;
            self.points[p].y += self.origin_point.y;
        }
        self.set_origin_point();
        self.set_normals();
    }



    pub fn set_normals(&mut self) {
        for i in 0..self.triangles.len() {
            self.triangles[i].normal = calc_normal([self.points[self.triangles[i].point1], self.points[self.triangles[i].point2], self.points[self.triangles[i].point3]]);
        }
    }
}

#[derive(Debug, Clone)]
pub struct Buffer {
    pub amt_of_points: usize,

    pub points: Vec<Point3>,
    pub triangles: Vec<Triangle>
}
impl Buffer {
    pub fn init_buffer() -> Buffer {
        Buffer {amt_of_points: 0, points: Vec::new(),triangles: Vec::new()}
    }

    pub fn load_mesh(&mut self, mesh: &Mesh) {
        for i in 0..mesh.points.len() {
            self.points.push(mesh.points[i]);
        }
        for i in 0..mesh.triangles.len() {
            let l = calc_normal([mesh.points[mesh.triangles[i].point1], mesh.points[mesh.triangles[i].point2], mesh.points[mesh.triangles[i].point3]]);
            self.triangles.push(Triangle { point1: mesh.triangles[i].point1+self.amt_of_points, point2: mesh.triangles[i].point2+self.amt_of_points, point3: mesh.triangles[i].point3+self.amt_of_points, origin_point: mesh.triangles[i].origin_point, color: mesh.triangles[i].color, normal: l });
        }
        self.amt_of_points += mesh.points.len();
    }
}

impl<'a> MeshAndBuffer for Buffer {
    fn sort_triangles_origin_point_z_distance_from_camera(&mut self, camera: &Camera, start: isize, end: isize) {
        if end <= start {return;}

        let mut i = start as usize;
        let mut j: isize = i as isize-1;
        let pivot = end as usize;

        while i < pivot {
            let length_point_i_to_camera = camera.origin_point.z - self.triangles[i].origin_point.z;
            let length_point_pivot_to_camera = camera.origin_point.z - self.triangles[pivot].origin_point.z;

            if length_point_i_to_camera > length_point_pivot_to_camera {
                j += 1;
                self.triangles.swap(i, j as usize);
            }
            i += 1;
        }
        j += 1;
        self.triangles.swap(j as usize, pivot);

        self.sort_triangles_origin_point_z_distance_from_camera(camera, start, j -1);
        self.sort_triangles_origin_point_z_distance_from_camera(camera, j + 1, end);
    }
    fn sort_triangles_in_mesh_origin_point_z_distance_from_camera(&mut self, camera: &Camera) {
        if self.triangles.len() == 0 {return;}

        let a = self.triangles.len()-1;
        self.sort_triangles_origin_point_z_distance_from_camera(camera, 0, a as isize);
    }

    fn sort_points_z_distance_from_camera(&mut self, camera: &Camera, start: isize, end: isize) {
        if end <= start {return;}

        let mut i = start as usize;
        let mut j: isize = i as isize-1;
        let pivot = end as usize;

        while i < pivot {
            let length_point_i_to_camera = camera.origin_point.z - self.points[i].z;
            let length_point_pivot_to_camera = camera.origin_point.z - self.points[pivot].z;

            if length_point_i_to_camera > length_point_pivot_to_camera {
                j += 1;
                self.points.swap(i, j as usize);
            }
            i += 1;
        }
        j += 1;
        self.points.swap(j as usize, pivot);

        self.sort_points_z_distance_from_camera(camera, start, j -1);
        self.sort_points_z_distance_from_camera(camera, j + 1, end);
    }
    fn sort_points_in_mesh_z_distance_from_camera(&mut self, camera: &Camera) {
        let a = self.points.len()-1;
        self.sort_points_z_distance_from_camera(camera, 0, a as isize);
    }
}

#[derive(Debug, Clone)]
pub struct Camera {
    pub origin_point: Point3,

    pub width: u32,
    pub height: u32,

    fov: f32,
    fov_radians: f32,

    near_clipping_plane: f32,
    far_clipping_plane: f32,

    pub forward_vec: Point3,
    pub right_vec: Point3,
    pub up_vec: Point3,
    pub look_vec: Point3,
    pub target_vec: Point3,
}

impl Camera {
    pub fn init(origin_point: Point3, width: u32, height: u32, fov: f32, near_clipping_plane: f32, far_clipping_plane: f32) -> Camera {
        Camera { 
            origin_point: (origin_point),
            width: (width),
            height: (height),
            fov: (fov),
            fov_radians: (fov*PI as f32/180.0),
            near_clipping_plane: (near_clipping_plane),
            far_clipping_plane: (far_clipping_plane),
            forward_vec: normalize_vector(Point3 { x: 0.0, y: 0.0, z: 1.0 }),
            right_vec: normalize_vector(Point3 { x: -1.0, y: 0.0, z: 0.0 }),
            up_vec: normalize_vector(Point3 { x: 0.0, y: 1.0, z: 0.0 }),
            look_vec: normalize_vector(Point3 { x: 0.0, y: 0.0, z: 1.0 }),
            target_vec: normalize_vector(Point3 { x: 0.0, y: 0.0, z: 1.0 }),
        }
    }

    pub fn move_camera(&mut self, x: f64, y: f64, z: f64) {
        self.origin_point.x += x;
        self.origin_point.y += y;
        self.origin_point.z += z;
    }
}

pub struct ObjLoader(pub Vec<Mesh>);
impl ObjLoader {
    pub fn init() -> ObjLoader {
        ObjLoader(Vec::new())
    }

    pub fn load_obj_file(&mut self, path:&str) -> Result<(), Box<dyn Error>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut amt_of_points = 0;
        let mut counter = 0;

        for line in reader.lines() {
            let line_string = line.unwrap();
            let mut line_split = line_string.split(" ");
            let line_0 = line_split.nth(0).unwrap().to_string(); 
            let obj_len = self.0.len();

            if line_0 == "o" {
                if counter > 0 {
                    amt_of_points += self.0[obj_len-1].points.len();
                }
                self.0.push(Mesh::init(line_split.nth(0).unwrap().to_string(), &[], &[], Color::RGB(180, 180, 180)));
                counter += 1;
            }
            
            if line_0 == "v" {
                let x = line_split.nth(0).unwrap().to_string();
                let y = line_split.nth(0).unwrap().to_string();
                let z = line_split.nth(0).unwrap().to_string();
                self.0[obj_len-1].points.push(Point3::init(x.parse()?, y.parse()?, z.parse()?));
            }

            if line_0 == "f" {
                /* */
                let mut p_1: u32 = line_split.nth(0).unwrap().to_string().parse()?;
                let mut p_2: u32 = line_split.nth(0).unwrap().to_string().parse()?;
                let mut p_3: u32 = line_split.nth(0).unwrap().to_string().parse()?;

                p_1 -= 1 + amt_of_points as u32;
                p_2 -= 1 + amt_of_points as u32;
                p_3 -= 1 + amt_of_points as u32;

                let p = Point3::init(
                    (self.0[obj_len-1].points[p_1 as usize].x + self.0[obj_len-1].points[p_2 as usize].x + self.0[obj_len-1].points[p_2 as usize].x)/3.0,
                    (self.0[obj_len-1].points[p_1 as usize].y + self.0[obj_len-1].points[p_2 as usize].y + self.0[obj_len-1].points[p_2 as usize].y)/3.0,
                    (self.0[obj_len-1].points[p_1 as usize].z + self.0[obj_len-1].points[p_2 as usize].z + self.0[obj_len-1].points[p_2 as usize].z)/3.0
                );
                let l = self.0[obj_len-1].color;
                let n = calc_normal([self.0[obj_len-1].points[p_1 as usize], self.0[obj_len-1].points[p_2 as usize], self.0[obj_len-1].points[p_3 as usize]]);
                self.0[obj_len-1].triangles.push(Triangle { point1: (p_1 as usize), point2: (p_2 as usize), point3: (p_3 as usize), origin_point: (p), color: (l), normal: n });
            }
        }  
        Ok(())
    }
}

pub trait SdlWrapper {
    fn draw_lines_w(&mut self, buffer: &mut Buffer, camera: &Camera);
    
    fn draw_triangle(&mut self, points: [Point; 3]);
    fn draw_triangles(&mut self, buffer: &mut Buffer, camera: &Camera);

    fn check_backface_culling(points: [Point3; 3], normal: Point3, camera: &Camera) -> bool;
}

impl SdlWrapper for Canvas<Window> {
    fn draw_lines_w(&mut self, buffer: &mut Buffer, camera: &Camera) {
        buffer.sort_triangles_in_mesh_origin_point_z_distance_from_camera(camera);
        for i in 0..buffer.triangles.len() {
            if buffer.triangles[i].normal.x * (buffer.points[buffer.triangles[i].point1].x - camera.origin_point.x) +
               buffer.triangles[i].normal.y * (buffer.points[buffer.triangles[i].point1].y - camera.origin_point.y) +
               buffer.triangles[i].normal.z * (buffer.points[buffer.triangles[i].point1].z - camera.origin_point.z) > 0.0
            {
                continue;
            }

            if buffer.triangles[i].color == None {
                self.set_draw_color(Color::RGB(255, 255, 255));
            }
            else {
                self.set_draw_color(Color::RGB(255, 255, 255));
            }
            let p1 = Point3::project_point(buffer.points[buffer.triangles[i].point1], camera);
            let p2 = Point3::project_point(buffer.points[buffer.triangles[i].point2], camera);
            let p3 =  Point3::project_point(buffer.points[buffer.triangles[i].point3], camera);

            if p1 == None || p2 == None || p3 == None {
                continue;
            }

            self.draw_line(p1.unwrap(), p2.unwrap()).unwrap();
            self.draw_line(p1.unwrap(), p3.unwrap()).unwrap();
            self.draw_line(p2.unwrap(), p3.unwrap()).unwrap();
        }
    }

    fn draw_triangle(&mut self, points_: [Point; 3]) {
        let mut points = points_;

        if points[1].y < points[0].y {points.swap(0, 1)};
        if points[2].y < points[0].y {points.swap(2, 0)};
        if points[2].y < points[1].y {points.swap(2, 1)};

        let m_0_2;
        let m_0_1;
        let m_1_2;

        if points[2].x as f32-points[0].x as f32 != 0.0 {
            m_0_2 = (points[2].y as f32-points[0].y as f32)/(points[2].x as f32-points[0].x as f32);
        }
        else {
            m_0_2 = (points[2].y as f32-points[0].y as f32)/(0.1);
        }

        if points[1].x as f32-points[0].x as f32 != 0.0 {
            m_0_1 = (points[1].y as f32-points[0].y as f32)/(points[1].x as f32-points[0].x as f32);
        }
        else {
            m_0_1 = (points[1].y as f32-points[0].y as f32)/(0.1);
        }

        if points[2].x as f32-points[1].x as f32 != 0.0 {
            m_1_2 = (points[2].y as f32-points[1].y as f32)/(points[2].x as f32-points[1].x as f32);
        }
        else {
            m_1_2 = (points[2].y as f32-points[1].y as f32)/(0.1);
        }
        

        let b_0_2 = points[0].y as f32 - m_0_2 * points[0].x as f32;
        let b_0_1 = points[0].y as f32 - m_0_1 * points[0].x as f32;
        let b_1_2 = points[1].y as f32 - m_1_2 * points[1].x as f32;

        for y in points[0].y+1..=points[2].y {
            if y <= points[1].y {
                let left_x = -((b_0_1-y as f32)/m_0_1);
                let right_x = -((b_0_2-y as f32)/m_0_2);
                self.draw_line(Point::new(left_x as i32, y), Point::new(right_x as i32, y)).unwrap(); 
            }
            else {
                let left_x = -((b_1_2-y as f32)/m_1_2);
                let right_x = -((b_0_2-y as f32)/m_0_2);
                self.draw_line(Point::new(left_x as i32, y), Point::new(right_x as i32, y)).unwrap(); 
            }
        }
    }

    fn draw_triangles(&mut self, buffer: &mut Buffer, camera: &Camera) {
        buffer.sort_triangles_in_mesh_origin_point_z_distance_from_camera(camera);
        for t in 0..buffer.triangles.len() {
            if Canvas::check_backface_culling([buffer.points[buffer.triangles[t].point1], buffer.points[buffer.triangles[t].point2], buffer.points[buffer.triangles[t].point3]], buffer.triangles[t].normal, camera) {
                continue;
            }
            
            if buffer.triangles[t].color == None {
                let mut r = 180 as f64 + (buffer.triangles[t].normal.x + buffer.triangles[t].normal.y + buffer.triangles[t].normal.z);
                let mut g = 180 as f64 + (buffer.triangles[t].normal.x + buffer.triangles[t].normal.y + buffer.triangles[t].normal.z);
                let mut buffer = 180 as f64 + (buffer.triangles[t].normal.x + buffer.triangles[t].normal.y + buffer.triangles[t].normal.z);
                
                if f64::abs(r) > 255.0 {
                    r = 180.0;
                }
                if f64::abs(g) > 255.0 {
                    g = 180.0;
                }
                if f64::abs(buffer) > 255.0 {
                    buffer = 180.0;
                }
                self.set_draw_color(Color::RGB(f64::abs(r) as u8, f64::abs(g) as u8, f64::abs(buffer) as u8));
            }
            else {
                let mut r = buffer.triangles[t].color.unwrap().r as f64 * (buffer.triangles[t].normal.x + buffer.triangles[t].normal.y + buffer.triangles[t].normal.z);
                let mut g = buffer.triangles[t].color.unwrap().g as f64 * (buffer.triangles[t].normal.x + buffer.triangles[t].normal.y + buffer.triangles[t].normal.z);
                let mut b = buffer.triangles[t].color.unwrap().b as f64 * (buffer.triangles[t].normal.x + buffer.triangles[t].normal.y + buffer.triangles[t].normal.z);
                
                if f64::abs(r) > 255.0 {
                    r = 240.0;
                }
                if f64::abs(g) > 255.0 {
                    g = 240.0;
                }
                if f64::abs(b) > 255.0 {
                    b = 240.0;
                }
                
                self.set_draw_color(Color::RGB(f64::abs(r) as u8, f64::abs(g) as u8, f64::abs(b) as u8));
            }
            let p1 = Point3::project_point(buffer.points[buffer.triangles[t].point1], camera);
            let p2 = Point3::project_point(buffer.points[buffer.triangles[t].point2], camera);
            let p3 = Point3::project_point(buffer.points[buffer.triangles[t].point3], camera);

            if p1 == None || p2 == None || p3 == None {
                continue;
            }

            self.draw_triangle([p1.unwrap(), p2.unwrap(), p3.unwrap()]);
            //self.present();
            //::std::thread::sleep(Duration::new(0, 1_000_000_00u32)); // sloppy FPS limit
        }
    }

    fn check_backface_culling(points: [Point3; 3], normal: Point3, camera: &Camera) -> bool {
        normal.x * (points[0].x - camera.origin_point.x) +
        normal.y * (points[0].y - camera.origin_point.y) +
        normal.z * (points[0].z - camera.origin_point.z) > 0.0
    }
}
