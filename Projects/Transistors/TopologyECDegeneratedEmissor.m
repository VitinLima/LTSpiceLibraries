function [Ib, Ic, Ie, I1, I2, Vb, Vc, Ve, Gm, Rpi, Vpi, Stp, er] = TopologyECDegeneratedEmissor(R1, R2, Rc, Re, Vcc, Beta, Is)
  Vt = 26e-3;
  Ib = 1e-5;
  lnIsBeta = log(Is/Beta);
  
  Vb = R2*(Vcc-Ib*R1)/(R1+R2);
  Ve = Ib*Re*(Beta+1);
  Vbe = Vb - Ve;
  er = lnIsBeta + Vbe/Vt - log(Ib);
  
  Stp = 1;
  StpLim = 100;
  
  ep = 1e-10;
  
  if sign(er) == 1
    Ibmin = Ib;
    while sign(er) == 1
      Ib *= 2;
      Vb = R2*(Vcc-Ib*R1)/(R1+R2);
      Ve = Ib*Re*(Beta+1);
      Vbe = Vb - Ve;
      er = lnIsBeta + Vbe/Vt - log(Ib);
    endwhile
    Ibmax = Ib;
  else
    Ibmax = Ib;
    while sign(er) == -1
      Ib /= 2;
      Vb = R2*(Vcc-Ib*R1)/(R1+R2);
      Ve = Ib*Re*(Beta+1);
      Vbe = Vb - Ve;
      er = lnIsBeta + Vbe/Vt - log(Ib);
    endwhile
    Ibmin = Ib;
  endif
  
  while and(abs(er) > abs(ep), Stp < StpLim)
    if sign(er) == 1
      Ibmin = Ib;
    else
      Ibmax = Ib;
    endif
    
    Ib = (Ibmax + Ibmin)/2;
    
    Vb = R2*(Vcc-Ib*R1)/(R1+R2);
    Ve = Ib*Re*(Beta+1);
    Vbe = Vb - Ve;
    er = lnIsBeta + Vbe/Vt - log(Ib);
  
    Stp++;
  endwhile
  
  Ic = Ib*Beta;
  Ie = Ib + Ic;
  I1 = (Vcc - Vb)/R1;
  I2 = Vb/R2;
  Vc = Vcc - Ic*Rc;
  Gm = Ic/Vt;
  Vpi = Vbe;
  Rpi = Vbe/Ib;
  
% Ex. 5.10 R1 = 16k, R2 = 9k, Rc = 1k, Re = 100, Vcc = 2.5V,
% Beta = 100, Is = 5e-17A, Que variação sofre a corrende de coletor se
% R1 for 1% maior q seu valor nominal?
% TopologyECDegeneratedEmissor(9e3,16e3,1e3,100,2.5,100,5e-17)
end