function j = jacobianPaMapping(rangeU)
  J2xi    = 0.5 * ( rangeU(2) - rangeU(1) );  
  j       = J2xi;
end