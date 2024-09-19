function [Ib, Ic, Ie, I1, I2, Vb, Vc, Ve, Gm, Rpi, Vpi, Stp, er] = TopologyECDegeneratedEmissor3(R1, R2, Rc, Re, Vcc, Beta, Is)
  Vt = 26e-3;
  lnIs = log(Is);
  lnBeta = log(Beta);
  Vbe = 800e-3;
  
  a = R2*R1/(R2+R1) + Re*(Beta+1);
  b = R2*Vcc/(R2+R1);
  
  Ib = (b - Vbe)/a;
  
  Vben = Vt*(lnBeta + log(Ib) - lnIs);
  
  Stp = 1;
  StpLim = 100;
  
  ep = 1e-10;
  er = Vbe - Vben;
  
  while and(abs(er) > abs(ep), Stp < StpLim)
    Vbe = Vben;
    Ib = (b - Vbe)/a;
    Vben = Vt*(lnBeta + log(Ib) - lnIs);
    
    er = Vbe - Vben;
    
    Stp++;
  endwhile
  
  Ic = Ib*Beta;
  Vc = Vcc - Ic*Rc;
  Ve = Ib*Re*(Beta+1);
  Ie = Ve/Re;
  Vb = Ve + Vbe;
  I2 = Vb/R2;
  I1 = (Vcc - Vb)/R1;
  Gm = Ic/Vt;
  Vpi = Vbe;
  Rpi = Vbe/Ib;
  
% Ex. 5.10 R1 = 16k, R2 = 9k, Rc = 1k, Re = 100, Vcc = 2.5V,
% Beta = 100, Is = 5e-17A, Que variação sofre a corrende de coletor se
% R1 for 1% maior q seu valor nominal?
% TopologyECDegeneratedEmissor(9e3,16e3,1e3,100,2.5,100,5e-17)
end