use svg::Document;
use svg::node::element::Path;
use svg::node::element::path::Data;
use num::complex::Complex;
use std::ops::Mul;

const EPSILON: f64 = 0.001;

#[derive(Debug)]
struct Mat {
    a: Complex<f64>,
    b: Complex<f64>,
    c: Complex<f64>,
    d: Complex<f64>,
}

impl Mat {
    fn new(a: Complex<f64>, b: Complex<f64>, c: Complex<f64>, d: Complex<f64>) -> Self {
        Mat {
            a:a,
            b:b,
            c:c,
            d:d,
        }
    }

    fn id() -> Self {
        Mat {
            a: Complex::new(1.0,0.0),
            b: Complex::new(0.0,0.0),
            c: Complex::new(0.0,0.0),
            d: Complex::new(1.0,0.0),
        }
    }

    fn adj(&self) -> Self {
        Mat {
            a: self.d,
            b: -self.b,
            c: -self.c,
            d: self.a,
        }
    }

    fn mob(&self, z: Complex<f64>) -> Complex<f64> {
        (self.a * z + self.b) / (self.c * z + self.d)
    }

    fn fix(&self) -> Complex<f64> {
        // gives the attracting fixed point
        // z = az+b/cz+d, with big cz+d
        // cz^2 + (d-a) z - b = 0
        let a = self.a;
        let b = self.b;
        let c = self.c;
        let d = self.d;
        if c.norm_sqr() == 0.0 {
            if a.norm_sqr() > d.norm_sqr() {
                Complex::new(1.0 / 0.0, 0.0)
            } else {
                b / (d-a)
            }
        } else {
            let disc = (d - a) * (d - a) + 4.0 * b * c;
            println!("{:?}", disc);
            let sd = if (a + d).re > 0.0 {
                -disc.sqrt()
            } else {
                disc.sqrt()
            };
            (a - d - sd) / (2.0 * c)
        }
    }
}

impl<'a,'b> Mul<&'b Mat> for &'a Mat {
    type Output = Mat;
    fn mul(self, rhs: &'b Mat) -> Mat {
        let v = &rhs;
        Mat {
            a: self.a * v.a + self.b * v.c,
            b: self.a * v.b + self.b * v.d,
            c: self.c * v.a + self.d * v.c,
            d: self.c * v.b + self.d * v.d,
        }
    }
}

impl Mul<Mat> for Mat {
    type Output = Mat;
    fn mul(self, rhs: Mat) -> Mat {
        &self * &rhs
    }
}

impl Mul<&Mat> for Mat {
    type Output = Mat;
    fn mul(self, rhs: &Mat) -> Mat {
        &self * rhs
    }
}



fn grandma(ta: Complex<f64>, tb: Complex<f64>) -> Kleinian {
    let i = Complex::i();
    let disc = ta * ta * tb * tb - 4.0 * ta * ta - 4.0 * tb * tb;
    let tab = 0.5 * (ta * tb - disc.sqrt());
    let scale = (tab - 2.0) * tb / (tb * tab - 2.0 * ta + 2.0 * i * tab);

    let a = Mat::new(ta / 2.0, (ta * tab - 2.0 * tb + 4.0 * i) / ((2.0 * tab + 4.0) * scale),
        scale * (ta * tab - 2.0 * tb - 4.0 * i) / (2.0 * tab - 4.0), ta / 2.0);
    let b = Mat::new((tb - 2.0 * i) / 2.0, tb / 2.0,
        tb / 2.0, (tb + 2.0 * i) / 2.0);
    return Kleinian::new(a,b);
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
enum Letter {
    A,
    B,
    AI,
    BI,
}

const A: Letter = Letter::A;
const B: Letter = Letter::B;
const AI: Letter = Letter::AI;
const BI: Letter = Letter::BI;

impl Letter {
    fn inv(&self) -> Self {
        match self {
            &A => AI,
            &B => BI,
            &AI => A,
            &BI => B,
        }
    }
}

struct Bag<T> {
    a: T,
    b: T,
    ainv: T,
    binv: T,
}

impl<T> Bag<T> {
    fn new(a: T, b: T, ainv: T, binv: T) -> Self {
        Bag {
            a: a,
            b: b,
            ainv: ainv,
            binv: binv,
        }
    }

    fn at(&self, l: Letter) -> &T {
        match l {
            A => &self.a,
            B => &self.b,
            AI => &self.ainv,
            BI => &self.binv,
        }
    }
}

// #[derive(Debug)]
struct Kleinian {
    mats: Bag<Mat>,
    data: Option<Data>,
    last: Complex<f64>,
}

impl Kleinian {
    fn new(a: Mat, b: Mat) -> Kleinian {
        let (ainv, binv) = (a.adj(), b.adj());
        let bag = Bag::new(a, b, ainv, binv);
        Kleinian {
            mats: bag,
            data: None,
            last: Complex::new(1.0, 0.0),
        }
    }

    fn mat(&self, l: Letter) -> &Mat {
        self.mats.at(l)
    }

    fn endfix(&self, l: Letter) -> Complex<f64> {
        let one = Complex::new(1.0, 0.0);
        match l {
            A => (&self.mats.binv * &self.mats.ainv).mob(one), // BAba
            B => self.mats.binv.mob(one), // aBAb
            AI => one, // baBA
            BI => self.mats.ainv.mob(one), // AbaB
        }
    }

    fn line(&mut self, z: Complex<f64>) {
        let data = self.data.take();
        self.data = match data {
            Some(d) => Some(d.line_to((z.re, z.im))),
            None => Some(Data::new().move_to((z.re, z.im))),
        };
        self.last = z;
        // mem::replace(&mut self.data, self.data.line_to((z.re, z.im)));
    }
}

fn branch(level: i64, l: Letter, t: &Mat, g: &mut Kleinian) {

    let (l1, l2, l3) = match l {
        A => (B, A, BI),
        B => (AI, B, A),
        AI => (BI, AI, B),
        BI => (A, BI, AI),
    };
    let one = Complex::new(1.0, 0.0);

    let t = t * &g.mat(l);
    let z = t.mob(g.endfix(l));
    // println!("{:?}", l);
    // println!("{:?}", z);

    if level <= 0 || (g.last - z).norm_sqr() < EPSILON * EPSILON {
        // println!("{:?}", z);
        g.line(z);
        return;
    }

    branch(level - 1, l1, &t, g);
    branch(level - 1, l2, &t, g);
    branch(level - 1, l3, &t, g);
}

fn limitset(level: i64, g: &mut Kleinian) {
    let one = Complex::new(1.0, 0.0);
    let t = Mat::id();
    g.line(one);
    branch(level - 1, A, &t, g);
    branch(level - 1, BI, &t, g);
    branch(level - 1, AI, &t, g);
    branch(level - 1, B, &t, g);
}

fn main() {
    // println!("{:?}", Mat::id());
    let one = Complex::new(1.0,0.0);
    // let zero = Complex::new(0.0,0.0);
    // let ma = Mat::new(one,one,zero,one);
    // let mb = Mat::new(one,zero,one,one);
    // let mr = &ma * &mb;
    // let ma_inv = ma.adj();
    // let mb_inv = mb.adj();
    // let mr2 = &ma_inv * &mb_inv;
    // let mr_inv = mr.adj();
    // let mr2_inv = mr2.adj();
    // println!("{:?}", (&mr * &mr2 * &mr_inv * &mr2_inv).mob(zero));

    // let mut g = grandma(Complex::new(1.73205080757,1.0), Complex::new(2.0,0.0));
    let mut g = grandma(Complex::new(2.0, 0.0), Complex::new(2.0, 0.0));
    // println!("{:?}", a);
    // println!("{:?}", b);
    // println!("{:?}", &a * &b);
    // println!("{:?}", &a * &b * &a.adj() * &b.adj());

    // let v = &g.b * &g.a * &g.binv * &g.ainv;
    // println!("{:?}", v);
    // println!("{:?}", v.mob(one));
    


    limitset(50, &mut g);

    let path = Path::new()
        .set("fill", "none")
        .set("stroke", "black")
        .set("stroke-width", 0.001)
        .set("d", g.data.unwrap());
    
    let document = Document::new()
        .set("viewBox", (-1.2, -1.2, 2.4, 2.4))
        .add(path);

    svg::save("image.svg", &document).unwrap();
}
