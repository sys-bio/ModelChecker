//! Infix text for expressions, for messages.

use sbml_io::{Constant, Math, Op};

/// Infix text for `m`, with only the parentheses precedence needs.
pub fn infix(m: &Math) -> String {
    let mut s = String::new();
    write(m, 0, &mut s);
    s
}

// Precedence levels: higher binds tighter.
const P_OR: u8 = 1;
const P_AND: u8 = 2;
const P_CMP: u8 = 3;
const P_ADD: u8 = 4;
const P_MUL: u8 = 5;
const P_NEG: u8 = 6;
const P_POW: u8 = 7;

fn num(x: f64) -> String {
    if x.fract() == 0.0 && x.abs() < 1e15 {
        format!("{}", x as i64)
    } else {
        format!("{x}")
    }
}

fn write(m: &Math, outer: u8, s: &mut String) {
    match m {
        Math::Number { value, .. } => {
            let v = value.value();
            if v < 0.0 && outer > P_ADD {
                s.push('(');
                s.push_str(&num(v));
                s.push(')');
            } else {
                s.push_str(&num(v));
            }
        }
        Math::Ident(id) => s.push_str(id),
        Math::Symbol(sym) => s.push_str(&sym.name),
        Math::Constant(c) => s.push_str(match c {
            Constant::True => "true",
            Constant::False => "false",
            Constant::Pi => "pi",
            Constant::ExponentialE => "exponentiale",
            Constant::NotANumber => "NaN",
            Constant::Infinity => "inf",
        }),
        Math::Apply(op, args) => apply(*op, args, outer, s),
        Math::Call(name, args) => call(name, args, s),
        Math::CallSymbol(sym, args) => call(&sym.name, args, s),
        Math::Root { degree, arg } => match degree {
            None => call("sqrt", std::slice::from_ref(arg), s),
            Some(d) => call("root", &[(**d).clone(), (**arg).clone()], s),
        },
        Math::Log { base, arg } => match base {
            None => call("log10", std::slice::from_ref(arg), s),
            Some(b) => call("log", &[(**b).clone(), (**arg).clone()], s),
        },
        Math::Piecewise { pieces, otherwise } => {
            s.push_str("piecewise(");
            let mut first = true;
            for (v, c) in pieces {
                if !first {
                    s.push_str(", ");
                }
                first = false;
                write(v, 0, s);
                s.push_str(", ");
                write(c, 0, s);
            }
            if let Some(o) = otherwise {
                if !first {
                    s.push_str(", ");
                }
                write(o, 0, s);
            }
            s.push(')');
        }
        Math::Lambda { args, body } => {
            s.push_str("lambda(");
            s.push_str(&args.join(", "));
            s.push_str(", ");
            write(body, 0, s);
            s.push(')');
        }
        Math::Semantics { math, .. } => write(math, outer, s),
        Math::Raw(_) => s.push_str("<unknown MathML>"),
    }
}

fn call(name: &str, args: &[Math], s: &mut String) {
    s.push_str(name);
    s.push('(');
    for (i, a) in args.iter().enumerate() {
        if i > 0 {
            s.push_str(", ");
        }
        write(a, 0, s);
    }
    s.push(')');
}

fn infix_op(op: Op) -> Option<(&'static str, u8)> {
    Some(match op {
        Op::Plus => (" + ", P_ADD),
        Op::Minus => (" - ", P_ADD),
        Op::Times => ("*", P_MUL),
        Op::Divide => ("/", P_MUL),
        Op::Power => ("^", P_POW),
        Op::Eq => (" == ", P_CMP),
        Op::Neq => (" != ", P_CMP),
        Op::Gt => (" > ", P_CMP),
        Op::Lt => (" < ", P_CMP),
        Op::Geq => (" >= ", P_CMP),
        Op::Leq => (" <= ", P_CMP),
        Op::And => (" && ", P_AND),
        Op::Or => (" || ", P_OR),
        _ => return None,
    })
}

fn apply(op: Op, args: &[Math], outer: u8, s: &mut String) {
    if op == Op::Minus && args.len() == 1 {
        let paren = outer > P_NEG;
        if paren {
            s.push('(');
        }
        s.push('-');
        write(&args[0], P_NEG + 1, s);
        if paren {
            s.push(')');
        }
        return;
    }
    if op == Op::Not && args.len() == 1 {
        s.push('!');
        write(&args[0], P_POW + 1, s);
        return;
    }
    match infix_op(op) {
        Some((text, prec)) if args.len() >= 2 => {
            let paren = outer > prec;
            if paren {
                s.push('(');
            }
            for (i, a) in args.iter().enumerate() {
                if i > 0 {
                    s.push_str(text);
                }
                // Left-associative: the right operand of - and / needs
                // parentheses at equal precedence; ^ is the other way round.
                let need = match (op, i) {
                    (Op::Power, 0) => prec + 1,
                    (Op::Power, _) => prec,
                    (Op::Minus | Op::Divide, i) if i > 0 => prec + 1,
                    _ => prec,
                };
                write(a, need, s);
            }
            if paren {
                s.push(')');
            }
        }
        // Empty sums and products, single-argument plus and times.
        Some(_) if args.len() == 1 => write(&args[0], outer, s),
        Some(_) if args.is_empty() => s.push_str(if op == Op::Times { "1" } else { "0" }),
        _ => call(op.tag(), args, s),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn id(s: &str) -> Math {
        Math::ident(s)
    }

    #[test]
    fn precedence() {
        let m = Math::Apply(
            Op::Divide,
            vec![
                Math::Apply(Op::Times, vec![id("Vm"), id("S")]),
                Math::Apply(Op::Plus, vec![id("Km"), id("S")]),
            ],
        );
        assert_eq!(infix(&m), "Vm*S/(Km + S)");
        let m = Math::Apply(
            Op::Minus,
            vec![id("a"), Math::Apply(Op::Minus, vec![id("b"), id("c")])],
        );
        assert_eq!(infix(&m), "a - (b - c)");
    }
}
