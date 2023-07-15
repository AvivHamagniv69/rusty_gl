use std::error::Error;
use std::f32::consts::PI;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::str;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::{AtomicU8, Ordering::Relaxed};
use std::sync::Mutex;
use std::thread;

use std::time::Instant;

use sdl2::{pixels::Color, rect::Point, render::Canvas, video::Window};

pub enum LightingMode {
    PerPixel,
    PerTriangle,
}

pub struct RustyGl {
    pub lighting_mode: LightingMode,
    backface_culling: AtomicBool,
    background_color: (AtomicU8, AtomicU8, AtomicU8, AtomicU8),
}
impl RustyGl {
    pub const fn init() -> RustyGl {
        RustyGl {
            lighting_mode: LightingMode::PerTriangle,
            background_color: (
                AtomicU8::new(0),
                AtomicU8::new(0),
                AtomicU8::new(0),
                AtomicU8::new(255),
            ),
            backface_culling: AtomicBool::new(true),
        }
    }

    pub fn get_backface_culling(&self) -> bool {
        self.backface_culling.load(Relaxed)
    }

    pub fn set_backface_culling(&self, val: bool) {
        self.backface_culling.store(val, Relaxed)
    }
}

pub static RUSTY_GL: RustyGl = RustyGl::init();

#[derive(Debug, Default, PartialEq, Clone, Copy, PartialOrd)]
pub struct Point3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

/*enum Axis {
    X(f32),
    Y(f32),
    Z(f32),
}*/

impl Point3 {
    /*fn find_point_on_line(line: [Point3; 2], known_axis: Axis) -> Result<Point3, ()> {
        let v = Point3 {x: line[0].x - line[1].x, y: line[0].y-line[1].y, z: line[0].z-line[1].z};
        match known_axis {
            Axis::X(x) => if v.z != 0.0 || v.y != 0.0 {

            },
            Axis::Y(y) => if v.x != 0.0 || v.z != 0.0 {

            },
            Axis::Z(z) => if v.x != 0.0 || v.y != 0.0 {
            },
        };
    }*/

    fn rotate_x(point_: Point3, origin_point: &Point3, degrees: f32) -> Point3 {
        let a = degrees * PI / 180.0;
        let mut point = point_;

        point.y -= origin_point.y;
        point.z -= origin_point.z;

        let y = point.y;
        let z = point.z;

        point.y = y * a.cos() - z * a.sin();
        point.z = z * a.sin() + z * a.cos();

        point.y += origin_point.y;
        point.z += origin_point.z;

        point
    }

    fn rotate_y(point_: Point3, origin_point: &Point3, degrees: f32) -> Point3 {
        let a = degrees * PI / 180.0;
        let mut point = point_;

        point.x -= origin_point.x;
        point.z -= origin_point.z;

        let x = point.x;
        let z = point.z;

        point.x = x * a.cos() + z * a.sin();
        point.z = -(x * a.sin()) + z * a.cos();

        point.x += origin_point.x;
        point.z += origin_point.z;
        point
    }

    fn rotate_z(point_: Point3, origin_point: &Point3, degrees: f32) -> Point3 {
        let a = degrees * PI / 180.0;
        let mut point = point_;

        point.x -= origin_point.x;
        point.y -= origin_point.y;

        let x = point.x;
        let y = point.y;

        point.x = x * a.cos() - y * a.sin();
        point.y = x * a.sin() + y * a.cos();

        point.x -= origin_point.x;
        point.y -= origin_point.y;
        point
    }

    fn project(point: Self, camera: &Camera) -> Option<Point3> {
        let width = camera.width as f32;
        let height = camera.height as f32;

        let perspective = [
            1.0 / ((width / height) * f32::tan(camera.fov_radians / 2.0)),
            1.0 / f32::tan(camera.fov_radians / 2.0),
            -((camera.far_clipping_plane + camera.near_clipping_plane)
                / (camera.far_clipping_plane - camera.near_clipping_plane)),
            -((2.0 * camera.far_clipping_plane * camera.near_clipping_plane)
                / (camera.far_clipping_plane - camera.near_clipping_plane)),
        ];

        let w = -point.z;

        if point.z < camera.near_clipping_plane || point.z > camera.far_clipping_plane {
            return None;
        }

        Some(Point3 {
            x: point.x * perspective[0] / w * camera.fov_radians * width + width / 2.0,
            y: point.y * perspective[1] / w * camera.fov_radians * height + height / 2.0,
            z: w,
        })
    }

    fn calculate_distance(point1: &Self, point2: &Self) -> f32 {
        f32::sqrt(
            (point1.x - point2.x).powi(2)
                + (point1.y - point2.y).powi(2)
                + (point1.z - point2.z).powi(2),
        )
    }

    fn get_dot_product(point1: &Point3, point2: &Point3) -> f32 {
        point1.x * point2.x + point1.y * point2.y + point1.z + point2.z
    }

    fn normalize_vector(point: Point3) -> Point3 {
        let m = (point.x.powi(2) + point.y.powi(2) + point.z.powi(2)).sqrt();

        Point3 {
            x: (point.x / m),
            y: (point.y / m),
            z: (point.z / m),
        }
    }
}

struct Texture {

}

#[derive(Debug, Clone, Copy)]
pub struct Triangle {
    point1: usize,
    point2: usize,
    point3: usize,

    origin_point: Point3,
    normal: Point3,

    color: Color,
    color_baked: Option<Color>,
}

impl Triangle {
    fn is_out_of_screen(points: [Point3; 3], width: u16, height: u16) -> bool {
        ((points[0].x > width as f32 || points[0].x < 0.0)
            && (points[1].x > width as f32 || points[1].x < 0.0)
            && (points[2].x > width as f32 || points[2].x < 0.0))
            || ((points[0].y > height as f32 || points[0].y < 0.0)
                && (points[1].y > height as f32 || points[1].y < 0.0)
                && (points[2].y > height as f32 || points[2].y < 0.0))
    }

    fn project(points: [Point3; 3], camera: &Camera) -> Option<TriangleProjected> {
        let p1 = match Point3::project(points[0], camera) {
            Some(x) => x,
            None => return None,
        };
        let p2 = match Point3::project(points[1], camera) {
            Some(x) => x,
            None => return None,
        };
        let p3 = match Point3::project(points[2], camera) {
            Some(x) => x,
            None => return None,
        };

        if Triangle::is_out_of_screen([p1, p2, p3], camera.width, camera.height) {
            return None;
        }

        Some(TriangleProjected {
            point1: p1,
            point2: p2,
            point3: p3,
            origin_point: None,
            normal: None,
            color: None,
        })
    }

    fn calc_normal(points: [Point3; 3]) -> Point3 {
        let vu = Point3 {
            x: points[1].x - points[0].x,
            y: points[1].y - points[0].y,
            z: points[1].z - points[0].z,
        };
        let vv = Point3 {
            x: points[2].x - points[0].x,
            y: points[2].y - points[0].y,
            z: points[2].z - points[0].z,
        };

        let a = Point3 {
            x: (vu.y * vv.z) - (vu.z * vv.y),
            y: (vu.z * vv.x) - (vu.x * vv.z),
            z: (vu.x * vv.y) - (vu.y * vv.x),
        };
        Point3::normalize_vector(a)
    }
}

#[allow(unused)]
pub struct TriangleProjected {
    pub point1: Point3,
    pub point2: Point3,
    pub point3: Point3,
    origin_point: Option<Point3>,
    normal: Option<Point3>,
    color: Option<Color>,
}

pub trait MeshAndBuffer {
    fn sort_triangles_in_mesh_origin_point_z_distance_from_camera(
        &mut self,
        camera: &Camera,
    ) -> Result<(), ()>;
}

#[derive(Debug, Default, Clone)]
pub struct Mesh {
    pub name: String,

    pub points: Vec<Point3>,
    pub triangles: Vec<Triangle>,

    pub origin_point: Point3,
}

impl MeshAndBuffer for Mesh {
    fn sort_triangles_in_mesh_origin_point_z_distance_from_camera(
        &mut self,
        camera: &Camera,
    ) -> Result<(), ()> {
        if self.triangles.is_empty() {
            return Err(());
        }

        let check = self
            .triangles
            .iter_mut()
            .try_for_each(|a| -> Result<(), ()> {
                a.origin_point = Point3 {
                    x: a.origin_point.x - camera.origin_point.x,
                    y: a.origin_point.y - camera.origin_point.y,
                    z: a.origin_point.z - camera.origin_point.z,
                };
                Ok(())
            });
        if check == Err(()) {
            return Err(());
        }

        self.triangles
            .sort_by(|a, b| b.origin_point.z.partial_cmp(&a.origin_point.z).unwrap());

        self.triangles.iter_mut().for_each(|a| {
            a.origin_point = Point3 {
                x: a.origin_point.x - camera.origin_point.x,
                y: a.origin_point.y - camera.origin_point.y,
                z: a.origin_point.z - camera.origin_point.z,
            }
        });
        Ok(())
    }
}

impl Mesh {
    pub fn init(name: String, points: &[Point3], triangles: &[Triangle]) -> Mesh {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;

        points.iter().for_each(|a| {
            x += a.x;
            y += a.y;
            z += a.z;
        });

        let p = Point3 {
            x: (x / points.len() as f32),
            y: (y / points.len() as f32),
            z: (z / points.len() as f32),
        };
        Mesh {
            name: (name),
            points: (points.to_vec()),
            triangles: (triangles.to_vec()),
            origin_point: (p),
        }
    }

    pub fn update_origin_point(&mut self) {
        let mut p = Point3 {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };

        self.points.iter().for_each(|a| {
            p.x += a.x;
            p.y += a.y;
            p.z += a.z;
        });

        p.x /= self.points.len() as f32;
        p.y /= self.points.len() as f32;
        p.z /= self.points.len() as f32;
        self.origin_point = p;
    }

    pub fn rotate_x(&mut self, mut degrees: f32) {
        degrees *= PI / 180.0;

        self.points.iter_mut().for_each(|a| {
            a.y -= self.origin_point.y;
            a.z -= self.origin_point.z;

            let y = a.y;
            let z = a.z;

            a.y = y * degrees.cos() - z * degrees.sin();
            a.z = y * degrees.sin() + z * degrees.cos();

            a.y += self.origin_point.y;
            a.z += self.origin_point.z;
        });

        self.update_normals();
    }

    pub fn rotate_y(&mut self, mut degrees: f32) {
        degrees *= PI / 180.0;

        self.points.iter_mut().for_each(|a| {
            a.x -= self.origin_point.x;
            a.z -= self.origin_point.z;

            let x = a.x;
            let z = a.z;

            a.x = x * degrees.cos() + z * degrees.sin();
            a.z = -x * degrees.sin() + z * degrees.cos();

            a.x += self.origin_point.x;
            a.z += self.origin_point.z;
        });
        self.update_normals();
    }

    pub fn rotate_z(&mut self, mut degrees: f32) {
        degrees *= PI / 180.0;

        self.points.iter_mut().for_each(|a| {
            a.x -= self.origin_point.x;
            a.y -= self.origin_point.y;

            let x = a.x;
            let y = a.y;

            a.x = x * degrees.cos() - y * degrees.sin();
            a.y = x * degrees.sin() + y * degrees.cos();

            a.x += self.origin_point.x;
            a.y += self.origin_point.y;
        });
        self.update_normals();
    }

    pub fn move_mesh(&mut self, x: f32, y: f32, z: f32) {
        self.points.iter_mut().for_each(|a| {
            a.x += x;
            a.y += y;
            a.z += z;
        });
        self.update_normals()
    }

    pub fn update_normals(&mut self) {
        self.triangles.iter_mut().for_each(|a| {
            a.normal = Triangle::calc_normal([
                self.points[a.point1],
                self.points[a.point2],
                self.points[a.point3],
            ])
        });
    }

    pub fn set_color(&mut self, color: Color) {
        self.triangles.iter_mut().for_each(|a| a.color = color);
    }

    pub fn bake(&mut self, lights: &[Light]) {
        self.triangles.iter_mut().for_each(|a| {
            a.color_baked = Some(Light::apply_shading_multiple_lights(
                &a.color, &a.normal, lights,
            ))
        })
    }

    pub fn delete_bake(&mut self) {
        self.triangles.iter_mut().for_each(|a| a.color_baked = None);
    }
}

#[derive(Debug, Default, Clone)]
pub struct Buffer {
    amt_of_points: usize,

    points: Vec<Point3>,
    triangles: Vec<Triangle>,
}
impl Buffer {
    pub fn init() -> Buffer {
        Buffer {
            amt_of_points: 0,
            points: Vec::new(),
            triangles: Vec::new(),
        }
    }

    pub fn load_mesh(&mut self, mesh: &Mesh) {
        mesh.points.iter().for_each(|a| self.points.push(*a));
        mesh.triangles.iter().for_each(|a| {
            self.triangles.push(Triangle {
                point1: a.point1 + self.amt_of_points,
                point2: a.point2 + self.amt_of_points,
                point3: a.point3 + self.amt_of_points,
                origin_point: a.origin_point,
                normal: a.normal,
                color: a.color,
                color_baked: a.color_baked,
            })
        });

        self.amt_of_points += mesh.points.len();
    }

    pub fn load_meshes(&mut self, meshes: &[Mesh]) {
        meshes.iter().for_each(|a| self.load_mesh(a));
    }
}

impl MeshAndBuffer for Buffer {
    fn sort_triangles_in_mesh_origin_point_z_distance_from_camera(
        &mut self,
        camera: &Camera,
    ) -> Result<(), ()> {
        if self.triangles.is_empty() {
            return Err(());
        }

        self.triangles.iter_mut().for_each(|a| {
            a.origin_point = Point3 {
                x: a.origin_point.x - camera.origin_point.x,
                y: a.origin_point.y - camera.origin_point.y,
                z: a.origin_point.z - camera.origin_point.z,
            };
        });

        self.triangles
            .sort_by(|a, b| b.origin_point.z.partial_cmp(&a.origin_point.z).unwrap());

        self.triangles.iter_mut().for_each(|a| {
            a.origin_point = Point3 {
                x: a.origin_point.x - camera.origin_point.x,
                y: a.origin_point.y - camera.origin_point.y,
                z: a.origin_point.z - camera.origin_point.z,
            }
        });
        Ok(())
    }
}

pub trait Lighting {
    fn to_light(&self) -> Light;
}
impl Lighting for Point3 {
    fn to_light(&self) -> Light {
        Light::Point(PointLight {
            origin_point: Point3 {
                x: self.x,
                y: self.y,
                z: self.z,
            },
            brightness: 1.0,
            color: Color {
                r: 255,
                g: 255,
                b: 255,
                a: 255,
            },
            coefficients: (1.0, 1.0, 1.0),
        })
    }
}
impl Lighting for Camera {
    fn to_light(&self) -> Light {
        Light::Point(PointLight {
            origin_point: self.origin_point,
            brightness: 1.0,
            color: Color {
                r: 255,
                g: 255,
                b: 255,
                a: 255,
            },
            coefficients: (1.0, 1.0, 1.0),
        })
    }
}

pub trait ConvertColors {
    fn convert_255_1(color: &Color) -> (f32, f32, f32, f32) {
        (
            color.r as f32 / 255.0,
            color.g as f32 / 255.0,
            color.b as f32 / 255.0,
            color.a as f32 / 255.0,
        )
    }

    fn convert_1_255(color: &(f32, f32, f32, f32)) -> Color {
        Color {
            r: (color.0 * 255.0) as u8,
            g: (color.1 * 255.0) as u8,
            b: (color.2 * 255.0) as u8,
            a: (color.3 * 255.0) as u8,
        }
    }
}
impl ConvertColors for Color {}

#[derive(Debug, Clone, Copy)]
pub struct PointLight {
    // TODO
    pub coefficients: (f32, f32, f32),

    pub origin_point: Point3,
    pub brightness: f32,
    pub color: Color,
}
impl PointLight {
    fn apply_shading(color: &Color, normal: &Point3, light: &PointLight) -> Color {
        let d = Point3::calculate_distance(normal, &light.origin_point).powi(2);
        if light.brightness / 10.0 >= d.sqrt() {
            return *color;
        }

        let c = Point3::get_dot_product(normal, &light.origin_point);

        let converted_color = Color::convert_255_1(color);
        let converted_light = Color::convert_255_1(&light.color);

        let c = (
            converted_color.0 * converted_light.0 * (c / d) * light.brightness,
            converted_color.1 * converted_light.1 * (c / d) * light.brightness,
            converted_color.2 * converted_light.2 * (c / d) * light.brightness,
            converted_color.3,
        );
        Color::convert_1_255(&c)
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Light {
    Point(PointLight),
}

impl Light {
    fn apply_shading_multiple_lights(color: &Color, normal: &Point3, lights: &[Light]) -> Color {
        let mut c = (0.0, 0.0, 0.0, 255.0);
        lights.iter().for_each(|a| {
            let temp = Color::convert_255_1(&Light::apply_shading(color, normal, a));
            c.0 += temp.0;
            c.1 += temp.1;
            c.2 += temp.2;
        });
        Color::convert_1_255(&c)
    }

    fn apply_shading(color: &Color, normal: &Point3, light: &Light) -> Color {
        match light {
            Light::Point(x) => PointLight::apply_shading(color, normal, x),
        }
    }

    pub fn unwrap_point(&self) -> Option<PointLight> {
        match self {
            Light::Point(x) => Some(x.clone()),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Camera {
    pub origin_point: Point3,

    pub width: u16,
    pub height: u16,

    fov_radians: f32,

    near_clipping_plane: f32,
    far_clipping_plane: f32,

    rotations: (f32, f32, f32),
}

impl Camera {
    pub fn init(
        origin_point: Point3,
        width: u16,
        height: u16,
        fov: f32,
        near_clipping_plane: f32,
        far_clipping_plane: f32,
    ) -> Camera {
        Camera {
            origin_point: (origin_point),
            width: (width),
            height: (height),
            fov_radians: (fov * PI / 180.0),
            near_clipping_plane: (near_clipping_plane),
            far_clipping_plane: (far_clipping_plane),
            rotations: (0.0, 0.0, 0.0),
        }
    }

    pub fn rotate_x(&mut self, degrees: f32) {
        self.rotations.0 += degrees;
        if self.rotations.0 == 360.0 {
            self.rotations.0 = 0.0
        } else if self.rotations.0 > 360.0 {
            self.rotations.0 -= 360.0
        }
    }

    pub fn rotate_y(&mut self, degrees: f32) {
        self.rotations.1 += degrees;
        if self.rotations.1 == 360.0 {
            self.rotations.1 = 0.0
        } else if self.rotations.1 > 360.0 {
            self.rotations.1 -= 360.0
        }
    }

    pub fn rotate_z(&mut self, degrees: f32) {
        self.rotations.2 += degrees;
        if self.rotations.2 == 360.0 {
            self.rotations.2 = 0.0
        } else if self.rotations.2 > 360.0 {
            self.rotations.2 -= 360.0
        }
    }

    fn perform_transformations(mut point: Point3, camera: &Camera) -> Point3 {
        point = Point3::rotate_x(point, &camera.origin_point, camera.rotations.0);
        point = Point3::rotate_y(point, &camera.origin_point, camera.rotations.1);
        point = Point3::rotate_z(point, &camera.origin_point, camera.rotations.0);
        point.x -= camera.origin_point.x;
        point.y -= camera.origin_point.y;
        point.z -= camera.origin_point.z;
        point
    }
}

#[derive(Debug)]
struct MyError;
impl std::fmt::Display for MyError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "There is an error: {}", self)
    }
}
impl Error for MyError {}

pub struct MeshLoader(pub Vec<Mesh>);
impl MeshLoader {
    pub fn load_obj_file(&mut self, path: &str) -> Result<(), Box<dyn Error>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut normals: Vec<Point3> = Vec::new();
        let mut amt_of_points = 0;
        let mut counter = 0;

        for line in reader.lines() {
            let line_string = line?;
            let mut line_split = line_string.split(' ');
            let line_0 = match line_split.next() {
                Some(x) => x,
                None => return Err(Box::new(MyError)),
            }
            .to_string();
            let obj_len = self.0.len();

            if line_0 == "o" {
                if counter > 0 {
                    amt_of_points += self.0[obj_len - 1].points.len();
                }
                self.0
                    .push(Mesh::init(line_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string(), &[], &[]));
                counter += 1;
            } else if line_0 == "vn" {
                let x = line_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string();
                let y = line_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string();
                let z = line_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string();
                normals.push(Point3 {
                    x: x.parse()?,
                    y: y.parse()?,
                    z: z.parse()?,
                });
            } else if line_0 == "v" {
                let x = line_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string();
                let y = line_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string();
                let z = line_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string();
                self.0[obj_len - 1].points.push(Point3 {
                    x: x.parse()?,
                    y: y.parse()?,
                    z: z.parse()?,
                });
            } else if line_0 == "f" {
                let mut p_1: u32;
                let mut p_2: u32;
                let mut p_3: u32;
                let mut normal: Option<Point3> = None;
                if line_string.contains("//") {
                    let part_1 = line_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string();
                    let part_2 = line_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string();
                    let part_3 = line_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string();

                    let mut part_1_split = part_1.split("//");
                    let mut part_2_split = part_2.split("//");
                    let mut part_3_split = part_3.split("//");

                    p_1 = part_1_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string().parse()?;
                    let _v_1: usize = part_1_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string().parse()?;

                    p_2 = part_2_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string().parse()?;
                    let v_2: usize = part_2_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string().parse()?;

                    p_3 = part_3_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string().parse()?;
                    let _v_3: usize = part_3_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string().parse()?;

                    normal = Some(normals[v_2 - 1]);
                } else {
                    p_1 = line_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string().parse()?;
                    p_2 = line_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string().parse()?;
                    p_3 = line_split.next().expect("cant read file, file is not in the correct format or has been tampered with").to_string().parse()?;
                }
                p_1 -= 1 + amt_of_points as u32;
                p_2 -= 1 + amt_of_points as u32;
                p_3 -= 1 + amt_of_points as u32;

                let p = Point3 {
                    x: (self.0[obj_len - 1].points[p_1 as usize].x
                        + self.0[obj_len - 1].points[p_2 as usize].x
                        + self.0[obj_len - 1].points[p_2 as usize].x)
                        / 3.0,
                    y: (self.0[obj_len - 1].points[p_1 as usize].y
                        + self.0[obj_len - 1].points[p_2 as usize].y
                        + self.0[obj_len - 1].points[p_2 as usize].y)
                        / 3.0,
                    z: (self.0[obj_len - 1].points[p_1 as usize].z
                        + self.0[obj_len - 1].points[p_2 as usize].z
                        + self.0[obj_len - 1].points[p_2 as usize].z)
                        / 3.0,
                };
                let l = Color {
                    r: 180,
                    g: 180,
                    b: 180,
                    a: 255,
                };
                let n;

                if normal.is_none() {
                    n = Triangle::calc_normal([
                        self.0[obj_len - 1].points[p_1 as usize],
                        self.0[obj_len - 1].points[p_2 as usize],
                        self.0[obj_len - 1].points[p_3 as usize],
                    ]);
                } else {
                    n = normal.unwrap();
                }
                self.0[obj_len - 1].triangles.push(Triangle {
                    point1: (p_1 as usize),
                    point2: (p_2 as usize),
                    point3: (p_3 as usize),
                    origin_point: p,
                    color: (l),
                    normal: n,
                    color_baked: None,
                });
            }
        }
        Ok(())
    }
}

fn backface_cull(first_point: &Point3, normal: &Point3) -> bool {
    normal.x * (first_point.x) + normal.y * (first_point.y) + normal.z * (first_point.z) >= 0.0
}

#[allow(unused_variables, unused_mut)]
pub trait SdlRendering {
    fn draw_lines_w(&mut self, buffer: &mut Buffer, camera: &Camera) -> Result<(), ()>;
    fn draw_triangle_single_threaded(&mut self, mut points_projected: [Point3; 3], color: Color) {}
}
impl SdlRendering for Canvas<Window> {
    fn draw_lines_w(&mut self, buffer: &mut Buffer, camera: &Camera) -> Result<(), ()> {
        let check = buffer.sort_triangles_in_mesh_origin_point_z_distance_from_camera(camera);
        if check == Err(()) {
            return Err(());
        }
        self.set_draw_color(Color::RGB(0, 0, 0));

        let check = buffer.triangles.iter().try_for_each(|a| {
            let p_trans1 = Camera::perform_transformations(buffer.points[a.point1], camera);
            let p_trans2 = Camera::perform_transformations(buffer.points[a.point2], camera);
            let p_trans3 = Camera::perform_transformations(buffer.points[a.point3], camera);
            let normal = Triangle::calc_normal([p_trans1, p_trans2, p_trans3]);

            if backface_cull(&p_trans1, &normal) {
                return Ok(());
            }
            let p1 = match Point3::project(buffer.points[a.point1], camera) {
                Some(x) => x,
                None => return Err(()),
            };
            let p2 = match Point3::project(buffer.points[a.point2], camera) {
                Some(x) => x,
                None => return Err(()),
            };
            let p3 = match Point3::project(buffer.points[a.point3], camera) {
                Some(x) => x,
                None => return Err(()),
            };

            let c1 = self.draw_line(
                Point::new(p1.x as i32, p1.y as i32),
                Point::new(p2.x as i32, p2.y as i32),
            );
            let c2 = self.draw_line(
                Point::new(p1.x as i32, p1.y as i32),
                Point::new(p3.x as i32, p3.y as i32),
            );
            let c3 = self.draw_line(
                Point::new(p2.x as i32, p2.y as i32),
                Point::new(p3.x as i32, p3.y as i32),
            );

            if c1 != Ok(()) || c2 != Ok(()) || c3 != Ok(()) {
                return Err(());
            }

            Ok(())
        });
        if check == Err(()) {
            return Err(());
        }
        Ok(())
    }

    fn draw_triangle_single_threaded(&mut self, mut points_projected: [Point3; 3], color: Color) {
        if points_projected[1].y < points_projected[0].y {
            points_projected.swap(0, 1)
        };
        if points_projected[2].y < points_projected[0].y {
            points_projected.swap(2, 0)
        };
        if points_projected[2].y < points_projected[1].y {
            points_projected.swap(2, 1)
        };

        let m_0_2;
        let m_0_1;
        let m_1_2;

        if points_projected[2].x - points_projected[0].x != 0.0 {
            m_0_2 = (points_projected[2].y - points_projected[0].y)
                / (points_projected[2].x - points_projected[0].x);
        } else {
            m_0_2 = (points_projected[2].y - points_projected[0].y) / (0.1);
        }

        if points_projected[1].x - points_projected[0].x != 0.0 {
            m_0_1 = (points_projected[1].y - points_projected[0].y)
                / (points_projected[1].x - points_projected[0].x);
        } else {
            m_0_1 = (points_projected[1].y - points_projected[0].y) / (0.1);
        }

        if points_projected[2].x - points_projected[1].x != 0.0 {
            m_1_2 = (points_projected[2].y - points_projected[1].y)
                / (points_projected[2].x - points_projected[1].x);
        } else {
            m_1_2 = (points_projected[2].y - points_projected[1].y) / (0.1);
        }

        let b_0_2 = points_projected[0].y - m_0_2 * points_projected[0].x;
        let b_0_1 = points_projected[0].y - m_0_1 * points_projected[0].x;
        let b_1_2 = points_projected[1].y - m_1_2 * points_projected[1].x;

        for y in points_projected[0].y as i32 + 1..=points_projected[1].y as i32 {
            let left_x = -((b_0_1 - y as f32) / m_0_1);
            let right_x = -((b_0_2 - y as f32) / m_0_2);

            for x in f32::min(left_x, right_x) as i32..=f32::max(left_x, right_x) as i32 {
                self.draw_point(Point::new(x, y)).unwrap();
            }
        }
        for y in points_projected[1].y as i32 + 1..=points_projected[2].y as i32 {
            let left_x = -((b_1_2 - y as f32) / m_1_2);
            let right_x = -((b_0_2 - y as f32) / m_0_2);

            for x in f32::min(left_x, right_x) as i32..=f32::max(left_x, right_x) as i32 {
                self.draw_point(Point::new(x, y)).unwrap();
            }
        }
    }
}

pub trait SdlWrapper {
    fn draw_triangles_single_threaded(
        &mut self,
        buffer: Buffer,
        camera: &mut Camera,
        lights: &Vec<Light>,
    ) -> Result<(), ()>;

    fn draw_triangles_multi_threaded(
        &mut self,
        buffer: Buffer,
        camera: &mut Camera,
        lights: &Vec<Light>,
        amt_threads: usize
    ) -> Result<(), ()>;
}

impl SdlWrapper for Canvas<Window> {
    fn draw_triangles_multi_threaded(
            &mut self,
            mut buffer: Buffer,
            camera: &mut Camera,
            lights: &Vec<Light>,
            amt_threads: usize
        ) -> Result<(), ()> {
            let start: Instant = Instant::now();
            let check = buffer.sort_triangles_in_mesh_origin_point_z_distance_from_camera(camera);
            if check == Err(()) {
                return Err(());
            }

        let triangles_mutex = Mutex::new(
            buffer
                .triangles
                .chunks(buffer.triangles.len() / amt_threads as usize),
        );
        let tri_vec: Mutex<Vec<(Point3, Point3, Point3, Color)>> = Mutex::new(vec![]);
        let points = buffer.points;

        thread::scope(|scope| {
            for _ in 0..amt_threads {
                scope.spawn(|| {
                    let mut vec: Vec<(Point3, Point3, Point3, Color)> = vec![];
                    let mut lock;
                    while let Some(chunk) = {
                        lock = triangles_mutex.lock().unwrap();
                        lock.next()
                    } {
                        drop(lock);
                        chunk.iter().for_each(|tri| {
                            let p_trans1 =
                                Camera::perform_transformations(points[tri.point1], camera);
                            let p_trans2 =
                                Camera::perform_transformations(points[tri.point2], camera);
                            let p_trans3 =
                                Camera::perform_transformations(points[tri.point3], camera);
                            let normal = Triangle::calc_normal([p_trans1, p_trans2, p_trans3]);

                            if backface_cull(&p_trans1, &normal) {
                                return;
                            }
                            let triangle =
                                match Triangle::project([p_trans1, p_trans2, p_trans3], camera) {
                                    Some(x) => x,
                                    None => return,
                                };

                            if !lights.is_empty() && tri.color_baked.is_none() {
                                vec.push((
                                    triangle.point1,
                                    triangle.point2,
                                    triangle.point3,
                                    Light::apply_shading_multiple_lights(
                                        &tri.color, &tri.normal, lights,
                                    ),
                                ))
                            } else if tri.color_baked.is_some() {
                                vec.push((
                                    triangle.point1,
                                    triangle.point2,
                                    triangle.point3,
                                    tri.color_baked.unwrap(),
                                ))
                            } else {
                                vec.push((
                                    triangle.point1,
                                    triangle.point2,
                                    triangle.point3,
                                    tri.color,
                                ))
                            }
                        });
                    }
                    tri_vec.lock().unwrap().append(&mut vec);
                });
            }
        });

        tri_vec.lock().unwrap().iter().for_each(|a| {
            self.set_draw_color(a.3);
            self.draw_triangle_single_threaded([a.0, a.1, a.2], a.3)
        });
        let duration: std::time::Duration = start.elapsed();
        println!("draw_triangles: {:#?}", duration);
        Ok(())
    }

    fn draw_triangles_single_threaded(
        &mut self,
        mut buffer: Buffer,
        camera: &mut Camera,
        lights: &Vec<Light>,
    ) -> Result<(), ()> {
        let start: Instant = Instant::now();
        let check = buffer.sort_triangles_in_mesh_origin_point_z_distance_from_camera(camera);
        if check == Err(()) {
            return Err(());
        }

        buffer.triangles.iter().for_each(|tri| {
            let p_trans1 = Camera::perform_transformations(buffer.points[tri.point1], camera);
            let p_trans2 = Camera::perform_transformations(buffer.points[tri.point2], camera);
            let p_trans3 = Camera::perform_transformations(buffer.points[tri.point3], camera);
            let normal = Triangle::calc_normal([p_trans1, p_trans2, p_trans3]);

            if RUSTY_GL.get_backface_culling() == true && backface_cull(&p_trans1, &normal) {
                return;
            }

            let triangle = match Triangle::project([p_trans1, p_trans2, p_trans3], camera) {
                Some(x) => x,
                None => return,
            };

            if !lights.is_empty() && tri.color_baked.is_none() {
                let color = Light::apply_shading_multiple_lights(
                    &tri.color,
                    &tri.normal,
                    lights,
                ); 
                self.set_draw_color(color);
                self.draw_triangle_single_threaded([
                    triangle.point1,
                    triangle.point2,
                    triangle.point3,
                ], color);
            } else if tri.color_baked.is_some() {
                self.set_draw_color(tri.color_baked.unwrap());
                self.draw_triangle_single_threaded([
                    triangle.point1,
                    triangle.point2,
                    triangle.point3,
                ], tri.color_baked.unwrap());
            } else {
                self.set_draw_color(tri.color);
                self.draw_triangle_single_threaded([
                    triangle.point1,
                    triangle.point2,
                    triangle.point3,],
                    tri.color
                );
            }
        });
        Ok(())
    }
}
