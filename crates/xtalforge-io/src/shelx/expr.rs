//! SHELX free-variable resolution and expression parsing (Task 3.3.5).
//!
//! `10m + p` is SHELX's numeric free-variable encoding; one-letter variables
//! `a`..`z` also name the `FVAR` values in order for simple expressions.

use xtalforge_core::Error;

pub(super) fn resolve_token(token: &str, fvar: &[f64]) -> Result<f64, Error> {
    if let Ok(value) = token.parse::<f64>() {
        return Ok(resolve_field(value, fvar));
    }
    Expr::evaluate(token, fvar)
}

/// SHELX's `10m + p` encoding: `m = 1` is a fixed value `p`, `m > 1`
/// multiplies the free variable `m` by `p`, and the negative form uses
/// `p * (fv(m) - 1)`.
pub(super) fn resolve_field(value: f64, fvar: &[f64]) -> f64 {
    if !value.is_finite() {
        return value;
    }
    let mode = (value.abs() / 10.0).floor() as usize;
    if mode == 0 {
        return value;
    }
    let sign = if value < 0.0 { -1.0 } else { 1.0 };
    let p = ((value - sign * 10.0 * mode as f64) * 1e8).round() / 1e8;
    if mode == 1 {
        return p;
    }
    let factor = fvar.get(mode - 1).copied().unwrap_or(1.0);
    if value > 0.0 {
        p * factor
    } else {
        p * (factor - 1.0)
    }
}

/// A small arithmetic parser for free-variable expressions: `+ - * /`,
/// parentheses, numbers and the one-letter variables `a`..`z`, which name
/// the `FVAR` values in order.
struct Expr<'a> {
    text: &'a [u8],
    at: usize,
    fvar: &'a [f64],
}

impl<'a> Expr<'a> {
    fn evaluate(text: &'a str, fvar: &'a [f64]) -> Result<f64, Error> {
        let mut parser = Self {
            text: text.as_bytes(),
            at: 0,
            fvar,
        };
        let value = parser.sum()?;
        parser.skip();
        if parser.at != parser.text.len() {
            return Err(expression_error(text));
        }
        if !value.is_finite() {
            return Err(expression_error(text));
        }
        Ok(value)
    }

    fn peek(&self) -> Option<u8> {
        self.text.get(self.at).copied()
    }

    fn bump(&mut self) -> Option<u8> {
        let byte = self.peek();
        if byte.is_some() {
            self.at += 1;
        }
        byte
    }

    fn skip(&mut self) {
        while self
            .peek()
            .is_some_and(|byte| byte == b' ' || byte == b'\t')
        {
            self.at += 1;
        }
    }

    fn sum(&mut self) -> Result<f64, Error> {
        let mut value = self.product()?;
        loop {
            self.skip();
            match self.peek() {
                Some(b'+') => {
                    self.at += 1;
                    value += self.product()?;
                }
                Some(b'-') => {
                    self.at += 1;
                    value -= self.product()?;
                }
                _ => return Ok(value),
            }
        }
    }

    fn product(&mut self) -> Result<f64, Error> {
        let mut value = self.unary()?;
        loop {
            self.skip();
            match self.peek() {
                Some(b'*') => {
                    self.at += 1;
                    value *= self.unary()?;
                }
                Some(b'/') => {
                    self.at += 1;
                    let divisor = self.unary()?;
                    if divisor == 0.0 {
                        return Err(Error::InvalidInput(
                            "SHELX free-variable expression divides by zero".to_owned(),
                        ));
                    }
                    value /= divisor;
                }
                _ => return Ok(value),
            }
        }
    }

    fn unary(&mut self) -> Result<f64, Error> {
        self.skip();
        match self.peek() {
            Some(b'-') => {
                self.at += 1;
                Ok(-self.unary()?)
            }
            Some(b'+') => {
                self.at += 1;
                self.unary()
            }
            _ => self.atom(),
        }
    }

    fn atom(&mut self) -> Result<f64, Error> {
        self.skip();
        match self.peek() {
            Some(b'(') => {
                self.at += 1;
                let value = self.sum()?;
                self.skip();
                if self.bump() != Some(b')') {
                    return Err(Error::InvalidInput(
                        "SHELX free-variable expression has an unclosed parenthesis".to_owned(),
                    ));
                }
                Ok(value)
            }
            Some(byte) if byte.is_ascii_digit() || byte == b'.' => self.number(),
            Some(byte) if byte.is_ascii_alphabetic() => self.variable(),
            _ => Err(Error::InvalidInput(
                "unreadable SHELX free-variable expression".to_owned(),
            )),
        }
    }

    fn number(&mut self) -> Result<f64, Error> {
        let mut digits = String::new();
        while let Some(byte) = self.peek() {
            if byte.is_ascii_digit() || byte == b'.' {
                digits.push(char::from(byte));
                self.at += 1;
            } else {
                break;
            }
        }
        digits
            .parse::<f64>()
            .map_err(|_| Error::InvalidInput(format!("bad SHELX number: {digits:?}")))
    }

    fn variable(&mut self) -> Result<f64, Error> {
        let Some(byte) = self.bump() else {
            return Err(Error::InvalidInput(
                "missing SHELX free-variable name".to_owned(),
            ));
        };
        let index = usize::from(byte.to_ascii_lowercase().saturating_sub(b'a'));
        self.fvar.get(index).copied().ok_or_else(|| {
            Error::InvalidInput(format!(
                "SHELX free variable {} has no FVAR value",
                char::from(byte)
            ))
        })
    }
}

fn expression_error(text: &str) -> Error {
    Error::InvalidInput(format!(
        "unreadable SHELX free-variable expression: {text:?}"
    ))
}
