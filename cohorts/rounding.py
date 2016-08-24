# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from decimal import Decimal, ROUND_HALF_UP

def round_float(f, digits, rounding=ROUND_HALF_UP):
    """
    Accurate float rounding from http://stackoverflow.com/a/15398691.
    """
    return Decimal(str(f)).quantize(Decimal(10) ** (-1 * digits),
                                    rounding=rounding)

def float_str(f, min_digits=2, max_digits=6):
    """
    Returns a string representing a float, where the number of
    significant digits is min_digits unless it takes more digits
    to hit a non-zero digit (and the number is 0 < x < 1).
    We stop looking for a non-zero digit after max_digits.
    """
    if f >= 1 or f <= 0:
        return str(round_float(f, min_digits))
    start_str = str(round_float(f, max_digits))
    digits = start_str.split(".")[1]
    non_zero_indices = []
    for i, digit in enumerate(digits):
        if digit != "0":
            non_zero_indices.append(i + 1)
    # Only saw 0s.
    if len(non_zero_indices) == 0:
        num_digits = min_digits
    else:
        # Of the non-zero digits, pick the num_digit'th of those (including any zeros)
        min_non_zero_indices = range(non_zero_indices[0], non_zero_indices[-1] + 1)[:min_digits]
        num_digits = min_non_zero_indices[-1]
    return str(round_float(f, num_digits))
