# This is the Polynomial Engine

***

 Specify a plynomial by its coefficients:
 
```p = Polynomial(1, 2, 5, name='g')```  
```print(p)```  

will show:  
```g(x) = 5x^2 + 2x + 1```

We can add and subtract polynomials by + and - respectively and 
```p.zeros``` gives us the complex roots of the polynomial.

```p.diff(n)```  
gives us the nth derivative.