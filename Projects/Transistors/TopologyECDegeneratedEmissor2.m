function [Ib, Ic, Ie, I1, I2, Vb, Vc, Ve, Gm, Rpi, Vpi, Stp, er] = TopologyECDegeneratedEmissor2(R1, R2, Rc, Re, Vcc, Beta, Is)
  Vt = 26e-3;
  
  Vb = R2*Vcc/(R1+R2);
  Vbe = 800e-3;
  Ve = Vb - Vbe;
  Ie = Ve/Re;
  Ic = Ie;
  
  Vbef = Vt*log(Ic/Is)
  er = Vbe - Vbef
  
  Stp = 1;
  StpLim = 100;
  
  ep = 1e-10;
  
  while and(abs(er) > abs(ep), Stp < StpLim)
    Vbe = Vbef;
    
    Ve = Vb - Vbe;
    Ie = Ve/Re;
    Ic = Ie;
    
    Vbef = Vt*log(Ic/Is);
    er = Vbe - Vbef
  
    Stp++
  endwhile
  
  Ib = Ic/Beta;
  I1 = Vcc/(R1+R2);
  I2 = I1;
  Vc = Vcc - Ic*Rc;
  Gm = Ic/Vt;
  Vpi = Vbe;
  Rpi = Vbe/Ib;
  
% Ex. 5.10 R1 = 16k, R2 = 9k, Rc = 1k, Re = 100, Vcc = 2.5V,
% Beta = 100, Is = 5e-17A, Que variação sofre a corrende de coletor se
% R1 for 1% maior q seu valor nominal?
% TopologyECDegeneratedEmissor(9e3,16e3,1e3,100,2.5,100,5e-17)
end