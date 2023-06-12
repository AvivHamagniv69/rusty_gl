use rand::Rng;
use sdl2::render::Canvas;
use std::f64::consts::PI;
use std::str;
use std::fs::File;
use std::error::Error;
use std::io::BufReader;
use std::io::BufRead;

use sdl2;
use sdl2::pixels::Color;
use sdl2::rect::Point;
use sdl2::video::Window;

#[derive(Debug, Clone, Copy)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,

    pub color: Option<Color>,
}
impl Point3 {
    pub fn init(x: f64, y: f64, z: f64, color: Option<Color>) -> Point3 {
        Point3 { x: (x), y: (y), z: (z), color: (color) }
    }    

    fn project_point(point: Point3, camera: &Camera) -> Point {
        let mut x: f64 = point.x-camera.origin_point.x;
        let mut y: f64 = point.y-camera.origin_point.y;
        let z: f64 = point.z-camera.origin_point.z;

        let t: f64 = 
        if camera.origin_point.z - z != 0.0 && camera.origin_point.z - z != -0.0 {
            (camera.near_clipping_plane as f64 - z) / (camera.origin_point.z - z)
        }
        else if camera.origin_point.z - z == -0.0 {
            (camera.near_clipping_plane as f64 - z) / (-0.000000000000000001)
        }
        else {
            (camera.near_clipping_plane as f64 - z) / (0.000000000000000001)
        };
        x = x + t * (camera.origin_point.x - x);
        y = y + t * (camera.origin_point.y - y);

        x *= camera.fov_radians as f64;
        y *= camera.fov_radians as f64;

        x *= camera.width as f64;
        y *= camera.height as f64;

        x += camera.width as f64/2.0;
        y += camera.height as f64/2.0;

        Point::new(x as i32, y as i32)
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

pub fn calc_origin_point_triangle(points: [Point3; 3]) -> Point3 {
    Point3 { 
        x: (points[0].x+points[1].x+points[2].x)/3.0,
        y: (points[0].y+points[1].y+points[2].y)/3.0,
        z: (points[0].z+points[1].z+points[2].z)/3.0,
        color: None 
    }
}

pub fn normalize_vector(point: Point3) -> Point3 {
    let m = f64::sqrt(f64::powi(point.x, 2) + f64::powi(point.y, 2)+ f64::powi(point.z, 2));

    Point3 { 
        x: (point.x/m),
        y: (point.y/m),
        z: (point.z/m),
        color: None 
    }
}

pub fn calc_normal(points: [Point3; 3]) -> Point3 {
    let vu = Point3::init(points[1].x-points[0].x, points[1].y-points[0].y, points[1].z-points[0].z, None);
    let vv = Point3::init(points[2].x-points[0].x, points[2].y-points[0].y, points[2].z-points[0].z, None);

    let a = Point3 {
        x: (vu.y*vv.z) - (vu.z*vv.y),
        y: (vu.z*vv.x) - (vu.x*vv.z),
        z: (vu.x*vv.y) - (vu.y*vv.x),
        color: None,
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

    origin_point: Option<Point3>,
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
        let p = Point3 {x: (x/points.len() as f64), y: (y/points.len() as f64), z: (z/points.len() as f64), color: (None)};
        Mesh {name: (name), points: (points.to_vec()), triangles: (triangles.to_vec()), origin_point: (Some(p)), color: (Some(color))}
    }

    pub fn set_origin_point(&mut self) {
        let mut p = Point3::init(0.0, 0.0, 0.0, None);
        for i in 0..self.points.len() {
            p.x += self.points[i].x;
            p.y += self.points[i].y;
            p.z += self.points[i].z;
        }
        p.x /= self.points.len() as f64;
        p.y /= self.points.len() as f64;
        p.z /= self.points.len() as f64;
        self.origin_point = Some(p);
    }

    pub fn rotate_y(&mut self, mut angle: f64) {
        angle *= PI/180.0;
        self.set_origin_point();
        for p in 0..self.points.len() {
            self.points[p].x = self.points[p].x * f64::cos(angle) + self.points[p].z * f64::sin(angle);
            self.points[p].z = -(self.points[p].x * f64:: sin(angle)) + self.points[p].z * f64::cos(angle);
        }
        self.set_origin_point();
    }

    pub fn rotate_z(&mut self, mut angle: f64) {
        angle *= PI/180.0;
        self.set_origin_point();
        for p in 0..self.points.len() {
            self.points[p].x = self.points[p].x * f64::cos(angle) - self.points[p].y * f64::sin(angle);
            self.points[p].y = self.points[p].x * f64:: sin(angle) + self.points[p].y * f64::cos(angle);
        }
        self.set_origin_point();
    }

    pub fn calc_normals(&mut self) {
        for i in 0..self.triangles.len() {
            self.triangles[i].normal = calc_normal([self.points[self.triangles[i].point1], self.points[self.triangles[i].point2], self.points[self.triangles[i].point3]]);
        }
    }
}

#[derive(Debug, Clone)]
pub struct Buffer<'a> {
    pub amt_of_points: usize,

    pub points: Vec<Point3>,
    pub triangles: Vec<Triangle>,

    pub active_camera: &'a Camera,
}
impl<'a> Buffer<'a> {
    pub fn init_buffer(camera: &Camera) -> Buffer {
        Buffer {amt_of_points: 0, points: Vec::new(),triangles: Vec::new(), active_camera: camera}
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

impl<'a> MeshAndBuffer for Buffer<'a> {
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
}

impl Camera {
    pub fn init(origin_point: Point3, width: u32, height: u32, fov: f32, near_clipping_plane: f32, far_clipping_plane: f32) -> Camera {
        Camera { origin_point: (origin_point), width: (width), height: (height), fov: (fov), fov_radians: (fov*PI as f32/180.0), near_clipping_plane: (near_clipping_plane), far_clipping_plane: (far_clipping_plane) }
    }
}

pub struct ObjLoader(pub Vec<Mesh>);
impl ObjLoader {
    pub fn init() -> ObjLoader {
        ObjLoader(Vec::new())
    }

    pub fn load_obj_file(&mut self, path:&str) -> Result<(), Box<dyn Error>> {
        let mut rng = rand::thread_rng();
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut amt_of_points = 0;

        for line in reader.lines() {
            let line_string = line.unwrap();
            let mut line_split = line_string.split(" ");
            let line_0 = line_split.nth(0).unwrap().to_string(); 
            let obj_len = self.0.len();

            if line_0 == "o" {
                if self.0.len() > 0 {
                    amt_of_points += self.0[obj_len-1].points.len();
                }
                self.0.push(Mesh::init(line_split.nth(0).unwrap().to_string(), &[], &[], Color::RGB(rng.gen_range(100..255), rng.gen_range(100..255), rng.gen_range(100..255))));
            }
            
            if line_0 == "v" {
                let x = line_split.nth(0).unwrap().to_string();
                let y = line_split.nth(0).unwrap().to_string();
                let z = line_split.nth(0).unwrap().to_string();
                let l = self.0[obj_len-1].color;
                self.0[obj_len-1].points.push(Point3::init(x.parse()?, y.parse()?, z.parse()?, l));
            }

            if line_0 == "f" {
                /* */
                let mut p_1: u32 = line_split.nth(0).unwrap().to_string().parse()?;
                let mut p_2: u32 = line_split.nth(0).unwrap().to_string().parse()?;
                let mut p_3: u32 = line_split.nth(0).unwrap().to_string().parse()?;

                p_1 -= 1 - amt_of_points as u32;
                p_2 -= 1 - amt_of_points as u32;
                p_3 -= 1 - amt_of_points as u32;

                let p = Point3::init(
                    (self.0[obj_len-1].points[p_1 as usize].x + self.0[obj_len-1].points[p_2 as usize].x + self.0[obj_len-1].points[p_2 as usize].x)/3.0,
                    (self.0[obj_len-1].points[p_1 as usize].y + self.0[obj_len-1].points[p_2 as usize].y + self.0[obj_len-1].points[p_2 as usize].y)/3.0,
                    (self.0[obj_len-1].points[p_1 as usize].z + self.0[obj_len-1].points[p_2 as usize].z + self.0[obj_len-1].points[p_2 as usize].z)/3.0,
                    None
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
    fn draw_points_w(&mut self, buffer: &Buffer, camera: &Camera);
    fn draw_lines_w(&mut self, buffer: &Buffer, camera: &Camera);
    
    fn draw_triangle(&mut self, points: [Point; 3]);
    fn draw_triangles(&mut self, buffer: &Buffer, camera: &Camera);

    fn draw_all(&mut self, buffer: &Buffer, camera: &Camera);
    fn check_backface_culling(points: [Point3; 3], normal: Point3, camera: &Camera) -> bool;
}

impl SdlWrapper for Canvas<Window> {
    fn draw_points_w(&mut self, buffer: &Buffer, camera: &Camera) {
        let mut b = buffer.clone();
        b.sort_points_in_mesh_z_distance_from_camera(camera);
        for p in b.points {
            if p.color == None {
                self.set_draw_color(Color::RGB(255, 255, 255));
            }
            else {
                self.set_draw_color(p.color.unwrap());
            }
            self.draw_point(Point3::project_point(p, camera)).unwrap();
        }
    }

    fn draw_lines_w(&mut self, buffer: &Buffer, camera: &Camera) {
        let mut b = buffer.clone();
        b.sort_triangles_in_mesh_origin_point_z_distance_from_camera(camera);
        for i in b.triangles {
            if i.normal.x * (b.points[i.point1].x - camera.origin_point.x) +
               i.normal.y * (b.points[i.point1].y - camera.origin_point.y) +
               i.normal.z * (b.points[i.point1].z - camera.origin_point.z) > 0.0
            {
                continue;
            }

            if i.color == None {
                self.set_draw_color(Color::RGB(255, 255, 255));
            }
            else {
                self.set_draw_color(Color::RGB(255, 255, 255));
            }
            let p1 = Point3::project_point(b.points[i.point1], camera);
            let p2 = Point3::project_point(b.points[i.point2], camera);
            let p3 =  Point3::project_point(b.points[i.point3], camera);
            self.draw_line(p1, p2).unwrap();
            self.draw_line(p1, p3).unwrap();
            self.draw_line(p2, p3).unwrap();
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

    fn draw_triangles(&mut self, buffer: &Buffer, camera: &Camera) {
        let mut b = buffer.clone();
        b.sort_triangles_in_mesh_origin_point_z_distance_from_camera(camera);
        for t in b.triangles {
            if Canvas::check_backface_culling([b.points[t.point1], b.points[t.point2], b.points[t.point3]], t.normal, camera) {
                continue;
            }
            
            if t.color == None {
                self.set_draw_color(Color::RGB(255, 255, 255));
            }
            else {
                self.set_draw_color(t.color.unwrap());
            }
            let p1 = Point3::project_point(buffer.points[t.point1], camera);
            let p2 = Point3::project_point(buffer.points[t.point2], camera);
            let p3 = Point3::project_point(buffer.points[t.point3], camera);
            self.draw_triangle([p1, p2, p3]);
        }
    }

    fn check_backface_culling(points: [Point3; 3], normal: Point3, camera: &Camera) -> bool {
        normal.x * (points[0].x - camera.origin_point.x) +
        normal.y * (points[0].y - camera.origin_point.y) +
        normal.z * (points[0].z - camera.origin_point.z) > 0.0
    }

    fn draw_all(&mut self, buffer: &Buffer, camera: &Camera) {
        
    }
}
