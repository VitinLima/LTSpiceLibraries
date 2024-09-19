function [Ib, Ic, Ie, I1, I2, Vb, Vc, Ve, Gm, Rpi, Vpi, Stp, er] = TopologyECDegeneratedEmissor4(R1, R2, Rc, Re, Vcc, Beta, Is)
  Vt = 26e-3;
  lnIs = log(Is);
  lnBeta = log(Beta);
  Vbe = 800e-3;
  Ib = 0;
  
  a = R2/(R2+R1);
  
  Vb = a*(Vcc - R1*Ib);
  
  Ve = Vb - Vbe;
  
  Ic = Ve/Re/(1+1/Beta);
  
  Ib = Ic/Beta;
  Vben = Vt*(log(Ic) - lnIs);
  
  Stp = 1;
  StpLim = 1e3;
  
  ep = 1e-10;
  er = Vbe - Vben;
  
##  ValVb = [Vb];
##  ValVe = [Ve];
##  ValVbe = [Vbe];
##  Valer = [er];
  
  while and(abs(er) > abs(ep), Stp < StpLim)
    
    Vbe = Vben;
    
    Vb = a*(Vcc - R1*Ib);
    
    Ve = Vb - Vbe;
    
    Ic = Ve/Re/(1+1/Beta);
    
    Ib = Ic/Beta;
    Vben = Vt*(log(Ic) - lnIs);
    
    er = Vbe - Vben;
  
##    ValVb(end+1) = Vb;
##    ValVe(end+1) = Ve;
##    ValVbe(end+1) = Vbe;
##    Valer(end+1) = er;
    
    Stp++;
  endwhile
  
  %[Ib, Ic, Ie, I1, I2, Vb, Vc, Ve, Gm, Rpi, Vpi, Stp, er]
  Vc = Vcc - Ic*Rc;
  
  Ie = Ic + Ib;
  I2 = Vb/R2;
  I1 = I2 + Ib;
  
  Gm = Ic/Vt;
  Vpi = Vbe;
  Rpi = Vbe/Ib;
  
##  figure;
##  hold on;
##  plot(ValVb, 'color', 'b');
##  plot(ValVe, 'color', 'g');
##  plot(ValVbe, 'color', 'r');
##  plot(Valer, 'color', 'b');
  
% Ex. 5.10 R1 = 16k, R2 = 9k, Rc = 1k, Re = 100, Vcc = 2.5V,
% Beta = 100, Is = 5e-17A, Que variação sofre a corrende de coletor se
% R1 for 1% maior q seu valor nominal?
% [Ib, Ic, Ie, I1, I2, Vb, Vc, Ve, Gm, Rpi, Vpi, Stp, er] = TopologyECDegeneratedEmissor4(16e3,9e3,1e3,100,2.5,100,5e-17)
end