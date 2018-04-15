
 Integer representation of floating point arithmetic suitable for FPGA designs
 
 Author: Yuri Gershtein 
 Date:   March 2018

 Functionality:

  *note* all integers are assumed to be signed

  all variables have units, stored in a map <string,int>, with string a unit (i.e. "phi") and int the power.
                   "2" is always present in the map, and it's int pair is referred to as 'shift'.
                   units are properly combined / propagated through calculations.
                   adding/subtracting variables with different units throws an exception.
                   adding/subtracting variables with different shifts is allowed and is handled correctly.

*Verilog_print()* method takes a vector of the outputs and produces a proper Verilog module

*calculate()* method re-calculates the variable double and int values based on its operands
                   returns false in case of overflows and/or mismatches between double and int calculations.

 the maximum and minimum values that the variable assumes are stored and updated each time calculate() is called
 if IMATH_ROOT is defined, all values are also stored in a histogram

*var_def*     (string name, string units, double fmax, double K):

                   define variable with bit value fval = K*ival, and maximum absolute value fmax.
                   calculates nbins on its own.
                   one can assign value to it using set_ methods.

*var_param*   (string name, double fval, int nbits):

                   define a parameter. K is calculated based on the fval and nbits.
		   one can assign value to it using set_ methods.

              (string name, string units, double fmax, double K):

                   define parameter with bit value fval = K*ival.
                   calculates nbins on its own.
                   one can assign value to it using set_ methods.

 *var_add*      (string name, var_base *p1, var_base *p2, double range = -1, int nmax = 18):
 
 *var_subtract* (string name, var_base *p1, var_base *p2, double range = -1, int nmax = 18):
 
                   add/subtract variables. Bit length increases by 1, but capped at nmax.
                   if range>0 specified, bit length is decreased to drop unnecessary high bits.

 *var_mult*    (string name, var_base *p1, var_base *p2, double range = -1, int nmax = 18):
 
                   multiplication. Bit length is a sum of the lengths of the operads, but capped at nmax.
                   if range>0 specified, bit length is decreased to drop unnecessary high bits or post-shift is reduced.

 *var_timesC*  (string name, var_base *p1, double cF, int ps = 17):
 
                   multiplication by a constant. Bit length stays the same
                   ps defines number of bits used to represent the constant

*var_DSP_postadd* (string name, var_base *p1, var_base *p2, var_base *p3, double range = -1, int nmax = 18):

                   explicit instantiation of the 3-clock DSP postaddition: p1*p2+p3
                   range and nmax have the same meaning as for the var_mult.

*var_shift*  (string name, var_base *p1, int shift):
 
                   shifts the variable right by shift (equivalent to multiplication by pow(2, -shift));
                   Units stay the same, nbits are adjusted.

 *var_neg*    (string name, var_base *p1):
 
                   multiplies the variable by -1

 *var_inv*    (string name, var_base *p1, double offset, int nbits, int n, unsigned int shift, mode m, int nbaddr=-1):
 
                   LUT-based inversion, f = 1./(offset + f1) and  i = 2^n / (offsetI + i1)
                   nbits is the width of the LUT (signed)
                   m is from enum mode {pos, neg, both} and refers to possible sign values of f
                            for pos and neg, the most significant bit of p1 (i.e. the sign bit) is ignored
                   shift is a shift applied in i1<->address conversions (used to reduce size of LUT)
                   nbaddr: if not specified, it is taken to be equal to p1->get_nbits()
                           

 *var_nounits* (string name, var_base *p1, int ps = 17):
 
                   convert a number with units to a number - needed for trig function expansion (i.e. 1 - 0.5*phi^2)
                   ps is a number of bits to represent the unit conversion constant

*var_adjustK* (string name, var_base *p1, double Knew, double epsilon = 1e-5, bool do_assert = false, int nbits = -1):

                   adjust variable shift so the K is as close to Knew as possible (needed for bit length adjustments) 
                   if do_assert is true, throw an exeption if Knew/Kold is not a power of two
                   epsilon is a comparison precision, nbits forces the bit length (possibly discarding MSBs)


*bool calculate* (int debug_level):

                     runs through the entire formula tree recalculating both ineteger and floating point values

                     returns true if everything is OK, false if obvious problems with the calculation exist, i.e
                                  -  integer value does not fit into the alotted number of bins
                                  -  integer value is more then 10% or more then 2 away from fval_/K_ 
                     debug_level:  0 - no warnings
                                   1 - limited warning
                                   2 - as 1, but also include explicit warnings when LUT was used out of its range
                                   3 - maximum complaints level


