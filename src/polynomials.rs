use std::collections::BTreeMap;

#[derive(Debug)]
struct Term {
    coefficient: i64,
    exponent: i64,
}

pub struct Polynomial {
    degree: u16,
    eqn: Vec<Term>,
}

impl Polynomial {
    pub fn print_polynomial(&self) {
        for term in self.eqn.iter() {
            print!("{}x^{} + ", term.coefficient, term.exponent)
        }
        println!()
    }
}

impl PartialEq for Term {
    fn eq(&self, other: &Self) -> bool {
        (self.coefficient == other.coefficient) && (self.exponent == other.exponent)
    }
}

pub fn multiply(poly1: &Polynomial, poly2: &Polynomial) -> Polynomial {
    let mut res_hash: BTreeMap<i64, i64> = BTreeMap::new();
    for term1 in poly1.eqn.iter() {
        for term2 in poly2.eqn.iter() {
            let res_exponent = term1.exponent + term2.exponent;
            let res_coefficient = term1.coefficient * term2.coefficient;

            res_hash
                .entry(res_exponent)
                .and_modify(|e| *e += res_coefficient)
                .or_insert(res_coefficient);

            // match res_hash.get(&res_exponent) {
            //     Some(cofficient) => res_hash.insert(res_exponent, cofficient + res_coefficient),
            //     None => res_hash.insert(res_exponent, res_coefficient),
            // };
        }
    }
    let mut res_eqn: Vec<Term> = Vec::new();
    for (key, value) in res_hash.iter().rev() {
        res_eqn.push(Term {
            coefficient: *value,
            exponent: *key,
        });
    }
    Polynomial {
        degree: poly1.degree + poly2.degree,
        eqn: res_eqn,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiply() {
        let poly1 = Polynomial {
            degree: 1,
            eqn: vec![
                Term {
                    coefficient: 4,
                    exponent: 1,
                },
                Term {
                    coefficient: -5,
                    exponent: 0,
                },
            ],
        };
        let poly2 = Polynomial {
            degree: 2,
            eqn: vec![
                Term {
                    coefficient: 2,
                    exponent: 2,
                },
                Term {
                    coefficient: 3,
                    exponent: 1,
                },
                Term {
                    coefficient: -6,
                    exponent: 0,
                },
            ],
        };
        // poly1.print_polynomial();
        // poly2.print_polynomial();
        let res = multiply(&poly1, &poly2);
        // res.print_polynomial();
        let expected_res_eqn = vec![
            Term {
                coefficient: 8,
                exponent: 3,
            },
            Term {
                coefficient: 2,
                exponent: 2,
            },
            Term {
                coefficient: -39,
                exponent: 1,
            },
            Term {
                coefficient: 30,
                exponent: 0,
            },
        ];
        assert_eq!(poly1.degree + poly2.degree, res.degree);
        assert_eq!(expected_res_eqn, res.eqn)
    }
}
