# maths
This project computes the smith normal form of a matrix over several principal ideal domains.

Supported domains:
- Integers (1)
- Gaussian integers (2)
- Polynomials with complex coefficients (3)
- Proper rational functions with complex coefficients (4)

The Smith Normal Form is useful in algebra for classifying finitely generated modules over a ring.  

Usage

Run:
python demo.py
Inputs: 
Matrices: Commas and semicolons.  Example for a 2x3 matrix 1 , 2 , 3 ; 4 , 5 , 6.  
Gaussian integer entries: use j.  Example 1+j with no space.  
Polynomial entries: Enter coefficients in ascending order with spaces to separate starting at x^0.  Example 1+x^2   is represented by 1 0 1.  Supports complex coefficients.  
Proper rational function entries: Enter numerator and denominator separated by /.  Example  ( 1+x^2 ) / ( 1-x^3 ) is represented by 1 0 1 / 1 0 0 -1.  Supports complex coefficients.  Please ensure entries are proper.  

