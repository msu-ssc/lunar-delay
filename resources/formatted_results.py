
from __future__ import annotations
 
import math
from dataclasses import dataclass
from decimal import Decimal

def format_result(value, uncertainty, unit="", sig_figs=None):
    """
    Round `value` and its 1-sigma `uncertainty` to a physically sensible
    number of digits
 
    sig_figs : how many significant figures to keep in the uncertainty.
        Default None -- auto-selects per the standard measurement-reporting
        convention: 1 significant figure, or 2 if the uncertainty's leading
        digit is 1 (dropping e.g. 1.4 to 1 throws away ~30% of the number;
        2 sig figs avoids that specifically for a leading 1). Pass an
        explicit int (e.g. sig_figs=2 or sig_figs=3) to always use that
        many significant figures instead of the auto rule.
 
    Either way, `value` is then rounded to that SAME decimal place as the
    (now-rounded) uncertainty. Never report a value more precisely than its
    own uncertainty supports.
 
    If `uncertainty` is zero, NaN, infinite, or negative-and-zero-after-abs
    (i.e. there's no usable error bar -- e.g. light_time_with_uncertainty()
    with n_samples=1), rounding-by-uncertainty is undefined, so `value` is
    returned unrounded (formatted to 6 significant figures as a fallback)
    and decimals/sig_figs are set to None so callers can detect this case.
 
    Uses Decimal(str(uncertainty)).adjusted() to find the uncertainty's
    order of magnitude -- exact, since it works off the same digit string
    Python would print, so there's no floating-point log10 edge case to
    guard against (the previous version needed an extra correction step
    for that; this doesn't).
    """
 
    digits = Decimal(str(abs(uncertainty)))
    exponent = digits.adjusted()  # order of magnitude of the leading digit
    if sig_figs is None:
        leading_digit = int(digits.as_tuple().digits[0])
        sig_figs = 2 if leading_digit == 1 else 1
 
    decimals = sig_figs - 1 - exponent  # digits after the decimal point; can be <= 0
    value=round(value, decimals)
    uncertainty=round(abs(uncertainty), decimals)
    return (value, uncertainty)
 
# a = format_result(1.325729024646947, 1.699267690872946e-07)
# print(a)