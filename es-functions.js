
/* SD from SD based on gains */
function sd_raw(sdg,r)
{
   return sdg/Math.sqrt(2*(1-r)) ;
}

/* Fisher's Zr transformation */
function Zr(r)
{
    return .5*Math.log((1+r)/(1-r)) ;
}

/* Inverse Fisher's Zr transformation */
function InvZr(z)
{
     return (Math.exp(2*z) - 1) / (Math.exp(2*z) + 1) ;
}

/* Compute d based on t-value and sample sizes */
function calc2()
{
  var frm = document.frm ;
  var t = parseFloat(frm.t.value) ;
  var grp1n = parseFloat(frm.grp1n.value) ;
  var grp2n = parseFloat(frm.grp2n.value) ;
  d = t*Math.sqrt((grp1n+grp2n)/(grp1n*grp2n));
  frm.d.value = round(d,4);
  v = vd(d,grp1n,grp2n);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d based on t-value and total sample sizes (assumes samples size
   are equal */
function calc3()
{
  var frm = document.frm ;
  var t      = parseFloat(frm.t.value) ;
  var totaln = parseFloat(frm.totaln.value) ;
  d = 2*t/Math.sqrt(totaln);
  frm.d.value = round(d,4);
  var grp1n = totaln*.5 ;
  var grp2n = totaln*.5 ;
  v = vd(d,grp1n,grp2n);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d based on a one-way F-test (two-groups) and sample sizes */
function calc4()
{
  var frm = document.frm ;
  var f = parseFloat(frm.f.value) ;
  var grp1n = parseFloat(frm.grp1n.value) ;
  var grp2n = parseFloat(frm.grp2n.value) ;
  d = Math.sqrt((f*(grp1n+grp2n))/(grp1n*grp2n));
  frm.d.value = round(d,4);
  v = vd(d,grp1n,grp2n);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d based on a one-way F-test (two-groups) and total sample sizes */
function calc5()
{
  var frm = document.frm ;
  var f = parseFloat(frm.f.value) ;
  var totaln = parseFloat(frm.totaln.value) ;
  var grp1n = totaln*.5 ;
  var grp2n = totaln*.5 ;
  d = 2*Math.sqrt((f/totaln));
  frm.d.value = round(d,4);
  v = vd(d,grp1n,grp2n);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d based on p-value of t and total sample size */
function calc6()
{
  var frm = document.frm ;
  var poft = parseFloat(frm.poft.value) ;
  var noft = parseFloat(frm.noft.value) ;
  df = noft - 2
  t = tfromp(df,poft) ;
  d = 2*t/Math.sqrt(noft);
  frm.d.value = round(d,4);
  var grp1n = noft*.5 ;
  var grp2n = noft*.5 ;
  v = vd(d,grp1n,grp2n);
  frm.v.value = round(v,6);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d based on p-value of t and separate sample sizes */
function calc7()
{
  var frm = document.frm ;
  var poft = parseFloat(frm.poft.value) ;
  var grp1n = parseFloat(frm.grp1n.value) ;
  var grp2n = parseFloat(frm.grp2n.value) ;
  df = grp1n+grp2n - 2 ;
  t = tfromp(df,poft) ;
  d = t*Math.sqrt((grp1n+grp2n)/(grp1n*grp2n));
  frm.d.value = round(d,4);
  v = vd(d,grp1n,grp2n);
  frm.v.value = round(v,6);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}


/* Compute d based on 2 by 2 frequency data */
function calc9()
{
  var frm = document.frm ;
  var grp1s = parseFloat(frm.grp1s.value) ;
  var grp1f = parseFloat(frm.grp1f.value) ;
  var grp2s = parseFloat(frm.grp2s.value) ;
  var grp2f = parseFloat(frm.grp2f.value) ;
  oddsratio = (grp1s*grp2f)/(grp2s*grp1f) ;
  loggedor  = Math.log(oddsratio) ;
  dl =  Math.sqrt(3)/Math.PI * loggedor ;
  dc = loggedor/1.65 ;
  vloggedor = 1/grp1s + 1/grp1f + 1/grp2s + 1/grp2f ;
  vl =  3/Math.pow(Math.PI,2) * vloggedor ;
  vc = vloggedor*.367 ;
  grp1ps = grp1s/(grp1s+grp1f) ;
  grp2ps = grp2s/(grp2s+grp2f) ;
  if (grp1ps<.5)  {     z1 = zfromp(grp1ps*2)*-1 ;    }
  else if (grp1ps>.5)  {     z1 = zfromp((1-grp1ps)/.5) ;    }
  else  {     z1=0 ;    }
  if (grp2ps<.5)     {     z2 = zfromp(grp2ps*2)*-1 ;    }
  else if (grp2ps>.5)    {     z2 = zfromp((1-grp2ps)/.5) ;    }
  else  {     z2=0 ;    }
  dp = z1 - z2 ;
  vp = (2*Math.PI*grp1ps*(1-grp1ps)* Math.exp(z1*z1)/(grp1s+grp1f))+(2*Math.PI*grp2ps*(1-grp2ps)* Math.exp(z2*z2)/(grp2s+grp2f));
  frm.dl.value =  round(dl,4);
  frm.dc.value =  round(dc,4);
  frm.dp.value =  round(dp,4);
  frm.vl.value = round(vl,6);
  frm.vc.value = round(vc,6);
  frm.vp.value = round(vp,6);
  frm.lowerl.value = round(lower_d(dl,vl),4) ;
  frm.upperl.value = round(upper_d(dl,vl),4) ;
  frm.lowerc.value = round(lower_d(dc,vc),4) ;
  frm.upperc.value = round(upper_d(dc,vc),4) ;
  frm.lowerp.value = round(lower_d(dp,vp),4) ;
  frm.upperp.value = round(upper_d(dp,vp),4) ;
}

/* Compute d based on 2 by 2 proportion data */
function calc10()
{
  var frm = document.frm ;
  var grp1ps = parseFloat(frm.grp1ps.value) ;
  var grp2ps = parseFloat(frm.grp2ps.value) ;
  var grp1n  = parseFloat(frm.grp1n.value) ;
  var grp2n  = parseFloat(frm.grp2n.value) ;
  oddsratio = (grp1ps*(1-grp2ps)) / (grp2ps*(1-grp1ps)) ;
  loggedor  = Math.log(oddsratio) ;
  dl =  Math.sqrt(3)/Math.PI * loggedor ;
  dc = loggedor/1.65 ;
  vloggedor = 1/(grp1ps*grp1n) + 1/((1-grp1ps)*grp1n) + 1/(grp2ps*grp2n) + 1/((1-grp2ps)*grp2n)
  vl =  3/Math.pow(Math.PI,2) * vloggedor ;
  vc = vloggedor*.367 ;
  if (grp1ps<.5)  {     z1 = zfromp(grp1ps*2)*-1 ;    }
  else if (grp1ps>.5)  {     z1 = zfromp((1-grp1ps)/.5) ;    }
  else  {     z1=0 ;    }
  if (grp2ps<.5)     {     z2 = zfromp(grp2ps*2)*-1 ;    }
  else if (grp2ps>.5)    {     z2 = zfromp((1-grp2ps)/.5) ;    }
  else  {     z2=0 ;    }
  dp = z1 - z2 ;
  vp = (2*Math.PI*grp1ps*(1-grp1ps)* Math.exp(z1*z1)/grp1n)+(2*Math.PI*grp2ps*(1-grp2ps)* Math.exp(z2*z2)/grp2n);
  frm.dl.value =  round(dl,4);
  frm.dc.value =  round(dc,4);
  frm.dp.value =  round(dp,4);
  frm.vl.value = round(vl,6);
  frm.vc.value = round(vc,6);
  frm.vp.value = round(vp,6);
  frm.lowerl.value = round(lower_d(dl,vl),4) ;
  frm.upperl.value = round(upper_d(dl,vl),4) ;
  frm.lowerc.value = round(lower_d(dc,vc),4) ;
  frm.upperc.value = round(upper_d(dc,vc),4) ;
  frm.lowerp.value = round(lower_d(dp,vp),4) ;
  frm.upperp.value = round(upper_d(dp,vp),4) ;
}

/* Compute d based on point-biserial r, unequal sample sizes */
function calc11()
{
  var frm = document.frm ;
  var r      = parseFloat(frm.r.value) ;
  var grp1n  = parseFloat(frm.grp1n.value) ;
  var grp2n  = parseFloat(frm.grp2n.value) ;
  p = grp1n/(grp1n+grp2n) ;
  d = r / Math.sqrt((1-r*r)*(p*(1-p))) ;
  v = vd(d,grp1n,grp2n) ;
  frm.d.value = round(d,4);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d based on point-biserial r, equal sample sizes */
function calc12()
{
  var frm = document.frm ;
  var r      = parseFloat(frm.r.value) ;
  var totaln = parseFloat(frm.totaln.value) ;
  var grp1n = totaln*.5 ;
  var grp2n = totaln*.5 ;
  d = 2*r / Math.sqrt(1-r*r) ;
  v = vd(d,grp1n,grp2n) ;
  frm.d.value = round(d,4);
  v = vd(d,grp1n,grp2n);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d based on phi coefficient from 2 by 2 */
function calc15()
{
  var frm = document.frm ;
  var r      = parseFloat(frm.r.value) ;
  var totaln = parseFloat(frm.totaln.value) ;
  d = 2*r / Math.sqrt(1-r*r) ;
  chisq = r*r * totaln ;
  v = (d*d)/chisq ;
  frm.d.value = round(d,4);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}


/* Compute d based on chi-square from 2 by 2 */
function calc17()
{
  var frm = document.frm ;
  var chisq  = parseFloat(frm.chisq.value) ;
  var totaln = parseFloat(frm.totaln.value) ;
  d = 2*Math.sqrt(chisq/(totaln-chisq)) ;
  v = (d*d)/chisq ;
  frm.d.value = round(d,4);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d based on chi-square from 2 by 2 */
function calc18()
{
  var frm = document.frm ;
  var pchisq  = parseFloat(frm.pchisq.value) ;
  var totaln = parseFloat(frm.totaln.value) ;
  chisq = Math.pow(zfromp(pchisq),2) ;
  d = 2*Math.sqrt(chisq/(totaln-chisq)) ;
  v = (d*d)/chisq ;
  frm.d.value = round(d,4) ;
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d from a frequency distribution */
function calc19()
{
  var frm = document.frm ;
  var grp1f = new Array() ;
  var grp2f = new Array() ;
  grp1f[0] = parseFloat(frm.grp1f0.value) ;
  grp2f[0] = parseFloat(frm.grp2f0.value) ;
  grp1f[1] = parseFloat(frm.grp1f1.value) ;
  grp2f[1] = parseFloat(frm.grp2f1.value) ;
  grp1f[2] = parseFloat(frm.grp1f2.value) ;
  grp2f[2] = parseFloat(frm.grp2f2.value) ;
  grp1f[3] = parseFloat(frm.grp1f3.value) ;
  grp2f[3] = parseFloat(frm.grp2f3.value) ;
  grp1f[4] = parseFloat(frm.grp1f4.value) ;
  grp2f[4] = parseFloat(frm.grp2f4.value) ;
  grp1f[5] = parseFloat(frm.grp1f5.value) ;
  grp2f[5] = parseFloat(frm.grp2f5.value) ;
  grp1f[6] = parseFloat(frm.grp1f6.value) ;
  grp2f[6] = parseFloat(frm.grp2f6.value) ;
  grp1f[7] = parseFloat(frm.grp1f7.value) ;
  grp2f[7] = parseFloat(frm.grp2f7.value) ;
  grp1f[8] = parseFloat(frm.grp1f8.value) ;
  grp2f[8] = parseFloat(frm.grp2f8.value) ;
  grp1f[9] = parseFloat(frm.grp1f9.value) ;
  grp2f[9] = parseFloat(frm.grp2f9.value) ;
  grp1n = 0 ;
  grp2n = 0 ;
  grp1fn = 0 ;
  grp2fn = 0 ;
  grp1f2n = 0 ;
  grp2f2n = 0 ;
  for (i=0;i<=9;i=i+1)
  {
   if (grp1f[i]>0)
     {
     grp1fn = grp1fn + i*grp1f[i] ;
     grp1n  = grp1n  + grp1f[i] ;
     grp1f2n = grp1f2n + i*i*grp1f[i] ;
     }
   if (grp2f[i]>0)
     {
     grp2fn = grp2fn + i*grp2f[i] ;
     grp2n  = grp2n  + grp2f[i] ;
     grp2f2n = grp2f2n + i*i*grp2f[i] ;
     }
  }
  grp1m = grp1fn/grp1n ;
  grp2m = grp2fn/grp2n ;
  grp1sd =   Math.sqrt(((grp1n * grp1f2n) - (grp1fn*grp1fn)) / (grp1n*grp1n)) ;
  grp2sd =   Math.sqrt(((grp2n * grp2f2n) - (grp2fn*grp2fn)) / (grp2n*grp2n)) ;
  sd_pooled = Math.sqrt((grp1sd*grp1sd*(grp1n-1)+grp2sd*grp2sd*(grp2n-1))/(grp1n+grp2n-2)) ;
  d = (grp1m-grp2m) / sd_pooled ;
  v = vd(d,grp1n,grp2n) ;
  frm.d.value = round(d,4);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d from a frequency distribution */
function calc20()
{
  var frm = document.frm ;
  var grp1n  = parseFloat(frm.grp1n.value) ;
  var grp2n  = parseFloat(frm.grp2n.value) ;
  var grp1f = new Array() ;
  var grp2f = new Array() ;
  var grp1p = new Array() ;
  var grp2p = new Array() ;
  grp1p[0] = parseFloat(frm.grp1p0.value) ;
  grp2p[0] = parseFloat(frm.grp2p0.value) ;
  grp1p[1] = parseFloat(frm.grp1p1.value) ;
  grp2p[1] = parseFloat(frm.grp2p1.value) ;
  grp1p[2] = parseFloat(frm.grp1p2.value) ;
  grp2p[2] = parseFloat(frm.grp2p2.value) ;
  grp1p[3] = parseFloat(frm.grp1p3.value) ;
  grp2p[3] = parseFloat(frm.grp2p3.value) ;
  grp1p[4] = parseFloat(frm.grp1p4.value) ;
  grp2p[4] = parseFloat(frm.grp2p4.value) ;
  grp1p[5] = parseFloat(frm.grp1p5.value) ;
  grp2p[5] = parseFloat(frm.grp2p5.value) ;
  grp1p[6] = parseFloat(frm.grp1p6.value) ;
  grp2p[6] = parseFloat(frm.grp2p6.value) ;
  grp1p[7] = parseFloat(frm.grp1p7.value) ;
  grp2p[7] = parseFloat(frm.grp2p7.value) ;
  grp1p[8] = parseFloat(frm.grp1p8.value) ;
  grp2p[8] = parseFloat(frm.grp2p8.value) ;
  grp1p[9] = parseFloat(frm.grp1p9.value) ;
  grp2p[9] = parseFloat(frm.grp2p9.value) ;
  grp1f[0] = grp1p[0]*grp1n ;
  grp2f[0] = grp2p[0]*grp2n ;
  grp1f[1] = grp1p[1]*grp1n ;
  grp2f[1] = grp2p[1]*grp2n ;
  grp1f[2] = grp1p[2]*grp1n ;
  grp2f[2] = grp2p[2]*grp2n ;
  grp1f[3] = grp1p[3]*grp1n ;
  grp2f[3] = grp2p[3]*grp2n ;
  grp1f[4] = grp1p[4]*grp1n ;
  grp2f[4] = grp2p[4]*grp2n ;
  grp1f[5] = grp1p[5]*grp1n ;
  grp2f[5] = grp2p[5]*grp2n ;
  grp1f[6] = grp1p[6]*grp1n ;
  grp2f[6] = grp2p[6]*grp2n ;
  grp1f[7] = grp1p[7]*grp1n ;
  grp2f[7] = grp2p[7]*grp2n ;
  grp1f[8] = grp1p[8]*grp1n ;
  grp2f[8] = grp2p[8]*grp2n ;
  grp1f[9] = grp1p[9]*grp1n ;
  grp2f[9] = grp2p[9]*grp2n ;
  grp1n = 0 ;
  grp2n = 0 ;
  grp1fn = 0 ;
  grp2fn = 0 ;
  grp1f2n = 0 ;
  grp2f2n = 0 ;
  for (i=0;i<=9;i=i+1)
  {
   if (grp1f[i]>0)
     {
     grp1fn = grp1fn + i*grp1f[i] ;
     grp1n  = grp1n  + grp1f[i] ;
     grp1f2n = grp1f2n + i*i*grp1f[i] ;
     }
   if (grp2f[i]>0)
     {
     grp2fn = grp2fn + i*grp2f[i] ;
     grp2n  = grp2n  + grp2f[i] ;
     grp2f2n = grp2f2n + i*i*grp2f[i] ;
     }
  }
  grp1m = grp1fn/grp1n ;
  grp2m = grp2fn/grp2n ;
  grp1sd =   Math.sqrt(((grp1n * grp1f2n) - (grp1fn*grp1fn)) / (grp1n*grp1n)) ;
  grp2sd =   Math.sqrt(((grp2n * grp2f2n) - (grp2fn*grp2fn)) / (grp2n*grp2n)) ;
  sd_pooled = Math.sqrt((grp1sd*grp1sd*(grp1n-1)+grp2sd*grp2sd*(grp2n-1))/(grp1n+grp2n-2)) ;
  d = (grp1m-grp2m) / sd_pooled ;
  v = vd(d,grp1n,grp2n) ;
  frm.d.value = round(d,4);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d based on means and full sample standard deviations */
function calc23()
{
  var frm = document.frm ;
  var grp1m = parseFloat(frm.grp1m.value) ;
  var grp1n = parseFloat(frm.grp1n.value) ;
  var grp2m = parseFloat(frm.grp2m.value) ;
  var grp2n = parseFloat(frm.grp2n.value) ;
  totaln = grp1n+grp2n ;
  var sdy = parseFloat(frm.sdy.value) ;
  sd_pooled = Math.sqrt(((sdy*sdy*(totaln-1)-((grp1m*grp1m+grp2m*grp2m-2*grp1m*grp2m)*(grp1n*grp2n))/totaln))/totaln) ;
  d = (grp1m-grp2m) / sd_pooled ;
  frm.d.value =  round(d,4);
  v = vd(d,grp1n,grp2n);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d based on gain score means */
function calc24()
{
  var frm = document.frm ;
  var grp1m = parseFloat(frm.grp1m.value) ;
  var grp2m = parseFloat(frm.grp2m.value) ;
  var grp1n = parseFloat(frm.grp1n.value) ;
  var grp2n = parseFloat(frm.grp2n.value) ;
  var grp1sdg = parseFloat(frm.grp1sdg.value) ;
  var grp2sdg = parseFloat(frm.grp2sdg.value) ;
  var grp1r = parseFloat(frm.grp1r.value) ;
  var grp2r = parseFloat(frm.grp2r.value) ;
  grp1sd = sd_raw(grp1sdg,grp1r) ;
  grp2sd = sd_raw(grp2sdg,grp2r) ;
  sd_pooled = Math.sqrt((grp1sd*grp1sd*(grp1n-1)+grp2sd*grp2sd*(grp2n-1))/(grp1n+grp2n-2)) ;
  d = (grp1m-grp2m) / sd_pooled ;
  v =  (2*(1-grp1r))/grp1n + (2*(1-grp2r))/grp2n + (d*d)/(2*(grp1n+grp2n)) ;
  frm.d.value =  round(d,4);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d based on gain score means */
function calc25()
{
  var frm = document.frm ;
  var grp1m = parseFloat(frm.grp1m.value) ;
  var grp2m = parseFloat(frm.grp2m.value) ;
  var grp1n = parseFloat(frm.grp1n.value) ;
  var grp2n = parseFloat(frm.grp2n.value) ;
  var grp1sd1 = parseFloat(frm.grp1sd1.value) ;
  var grp2sd1 = parseFloat(frm.grp2sd1.value) ;
  var grp1sd2 = parseFloat(frm.grp1sd2.value) ;
  var grp2sd2 = parseFloat(frm.grp2sd2.value) ;
  var grp1t = parseFloat(frm.grp1t.value) ;
  var grp2t = parseFloat(frm.grp2t.value) ;
  grp1sd = Math.sqrt(((grp1sd1*grp1sd1) + (grp1sd2*grp1sd2))/2) ;
  grp2sd = Math.sqrt(((grp2sd1*grp2sd1) + (grp2sd2*grp2sd2))/2) ;
  sd_pooled = Math.sqrt((grp1sd*grp1sd*(grp1n-1)+grp2sd*grp2sd*(grp2n-1))/(grp1n+grp2n-2)) ;
  d = (grp1m-grp2m) / sd_pooled ;
  grp1r = Math.abs(((grp1m*grp1m*grp1n)-(grp1sd1*grp1sd1*grp1t*grp1t + grp1sd2*grp1sd2*grp1t*grp1t))/(2*grp1sd1*grp1sd2*grp1t*grp1t)) ;
  grp2r = Math.abs(((grp2m*grp2m*grp2n)-(grp2sd1*grp2sd1*grp2t*grp2t + grp2sd2*grp2sd2*grp2t*grp2t))/(2*grp2sd1*grp2sd2*grp2t*grp2t)) ;
  v =  (2*(1-grp1r))/grp1n + (2*(1-grp2r))/grp2n + (d*d)/(2*(grp1n+grp2n)) ;
  frm.d.value =  round(d,4);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d based on gain score means */
function calc26()
{
  var frm = document.frm ;
  var grp1m = parseFloat(frm.grp1m.value) ;
  var grp2m = parseFloat(frm.grp2m.value) ;
  var grp1n = parseFloat(frm.grp1n.value) ;
  var grp2n = parseFloat(frm.grp2n.value) ;
  var grp1sd1 = parseFloat(frm.grp1sd1.value) ;
  var grp2sd1 = parseFloat(frm.grp2sd1.value) ;
  var grp1sd2 = parseFloat(frm.grp1sd2.value) ;
  var grp2sd2 = parseFloat(frm.grp2sd2.value) ;
  var grp1r = parseFloat(frm.grp1r.value) ;
  var grp2r = parseFloat(frm.grp2r.value) ;
  grp1sd = Math.sqrt(((grp1sd1*grp1sd1) + (grp1sd2*grp1sd2))/2) ;
  grp2sd = Math.sqrt(((grp2sd1*grp2sd1) + (grp2sd2*grp2sd2))/2) ;
  sd_pooled = Math.sqrt((grp1sd*grp1sd*(grp1n-1)+grp2sd*grp2sd*(grp2n-1))/(grp1n+grp2n-2)) ;
  d = (grp1m-grp2m) / sd_pooled ;
  v =  (2*(1-grp1r))/grp1n + (2*(1-grp2r))/grp2n + (d*d)/(2*(grp1n+grp2n)) ;
  frm.d.value =  round(d,4);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

function calc27()
{
  var frm = document.frm ;
  var grp1m  = new Array() ;
  var grp1sd = new Array() ;
  var grp1n  = new Array() ;
  var grp2m  = new Array() ;
  var grp2sd = new Array() ;
  var grp2n  = new Array() ;
  grp1m[0] = parseFloat(frm.grp1s0m.value) ;
  grp1m[1] = parseFloat(frm.grp1s1m.value) ;
  grp1m[2] = parseFloat(frm.grp1s2m.value) ;
  grp1m[3] = parseFloat(frm.grp1s3m.value) ;
  grp1m[4] = parseFloat(frm.grp1s4m.value) ;
  grp1m[5] = parseFloat(frm.grp1s5m.value) ;
  grp1m[6] = parseFloat(frm.grp1s6m.value) ;
  grp1m[7] = parseFloat(frm.grp1s7m.value) ;
  grp1m[8] = parseFloat(frm.grp1s8m.value) ;
  grp1sd[0] = parseFloat(frm.grp1s0sd.value) ;
  grp1sd[1] = parseFloat(frm.grp1s1sd.value) ;
  grp1sd[2] = parseFloat(frm.grp1s2sd.value) ;
  grp1sd[3] = parseFloat(frm.grp1s3sd.value) ;
  grp1sd[4] = parseFloat(frm.grp1s4sd.value) ;
  grp1sd[5] = parseFloat(frm.grp1s5sd.value) ;
  grp1sd[6] = parseFloat(frm.grp1s6sd.value) ;
  grp1sd[7] = parseFloat(frm.grp1s7sd.value) ;
  grp1sd[8] = parseFloat(frm.grp1s8sd.value) ;
  grp1n[0] = parseFloat(frm.grp1s0n.value) ;
  grp1n[1] = parseFloat(frm.grp1s1n.value) ;
  grp1n[2] = parseFloat(frm.grp1s2n.value) ;
  grp1n[3] = parseFloat(frm.grp1s3n.value) ;
  grp1n[4] = parseFloat(frm.grp1s4n.value) ;
  grp1n[5] = parseFloat(frm.grp1s5n.value) ;
  grp1n[6] = parseFloat(frm.grp1s6n.value) ;
  grp1n[7] = parseFloat(frm.grp1s7n.value) ;
  grp1n[8] = parseFloat(frm.grp1s8n.value) ;
  grp2m[0] = parseFloat(frm.grp2s0m.value) ;
  grp2m[1] = parseFloat(frm.grp2s1m.value) ;
  grp2m[2] = parseFloat(frm.grp2s2m.value) ;
  grp2m[3] = parseFloat(frm.grp2s3m.value) ;
  grp2m[4] = parseFloat(frm.grp2s4m.value) ;
  grp2m[5] = parseFloat(frm.grp2s5m.value) ;
  grp2m[6] = parseFloat(frm.grp2s6m.value) ;
  grp2m[7] = parseFloat(frm.grp2s7m.value) ;
  grp2m[8] = parseFloat(frm.grp2s8m.value) ;
  grp2sd[0] = parseFloat(frm.grp2s0sd.value) ;
  grp2sd[1] = parseFloat(frm.grp2s1sd.value) ;
  grp2sd[2] = parseFloat(frm.grp2s2sd.value) ;
  grp2sd[3] = parseFloat(frm.grp2s3sd.value) ;
  grp2sd[4] = parseFloat(frm.grp2s4sd.value) ;
  grp2sd[5] = parseFloat(frm.grp2s5sd.value) ;
  grp2sd[6] = parseFloat(frm.grp2s6sd.value) ;
  grp2sd[7] = parseFloat(frm.grp2s7sd.value) ;
  grp2sd[8] = parseFloat(frm.grp2s8sd.value) ;
  grp2n[0] = parseFloat(frm.grp2s0n.value) ;
  grp2n[1] = parseFloat(frm.grp2s1n.value) ;
  grp2n[2] = parseFloat(frm.grp2s2n.value) ;
  grp2n[3] = parseFloat(frm.grp2s3n.value) ;
  grp2n[4] = parseFloat(frm.grp2s4n.value) ;
  grp2n[5] = parseFloat(frm.grp2s5n.value) ;
  grp2n[6] = parseFloat(frm.grp2s6n.value) ;
  grp2n[7] = parseFloat(frm.grp2s7n.value) ;
  grp2n[8] = parseFloat(frm.grp2s8n.value) ;
  // create some needed variables
  grp1nm = 0 ;
  grp2nm = 0 ;
  grp1nm2 = 0 ;
  grp2nm2 = 0 ;
  grp1totaln = 0 ;
  grp2totaln = 0 ;
  grp1v = 0 ;
  grp2v = 0 ;
  k1 = 0 ;
  k2 = 0 ;
  // determine which radio button is checked
  if (frm.sdtype[0].checked)
  {
    var sdtype = 1 ;
  }
  else if (frm.sdtype[1].checked)
  {
    var sdtype = 2 ;
  }

  for (i=0;i<=8;i=i+1)
  {
   if (grp1n[i]>0 && grp1sd[i]>0 && grp1m[i]!="")
   {
     k1 = k1 + 1 ;
     grp1nm = grp1nm + grp1m[i]*grp1n[i] ;
     grp1nm2 = grp1nm2 + grp1m[i]*grp1m[i]*grp1n[i] ;
     grp1totaln = grp1totaln + grp1n[i] ;
     grp1v = grp1v  + (grp1sd[i]*grp1sd[i]*(grp1n[i]-1)) ;

   }
   if (grp2n[i]>0 && grp2sd[i]>0 && grp2m[i]!="")
   {
     k2 = k2 + 1 ;
     grp2nm = grp2nm + grp2m[i]*grp2n[i] ;
     grp2nm2 = grp2nm2 + grp2m[i]*grp2m[i]*grp2n[i] ;
     grp2totaln = grp2totaln + grp2n[i] ;
     grp2v = grp2v  + (grp2sd[i]*grp2sd[i]*(grp2n[i]-1)) ;
   }
  }

  grp1m = grp1nm/grp1totaln ;
  grp2m = grp2nm/grp2totaln ;
  grp1sd =  Math.sqrt(grp1v/(grp1totaln-k1)) ;
  grp2sd =  Math.sqrt(grp2v/(grp2totaln-k2)) ;
  if (sdtype==2)
  {
    ssb1 = (grp1nm2 - (grp1nm*grp1nm)/grp1totaln) ;
    ssb2 = (grp2nm2 - (grp2nm*grp2nm)/grp2totaln) ;
    grp1sd =  Math.sqrt((ssb1+grp1v)/(grp1totaln-1)) ;
    grp2sd =  Math.sqrt((ssb2+grp2v)/(grp2totaln-1)) ;
  }
  sd_pooled = Math.sqrt((grp1sd*grp1sd*(grp1totaln-1)+grp2sd*grp2sd*(grp2totaln-1))/(grp1totaln+grp2totaln-2)) ;
  d = (grp1m-grp2m) / sd_pooled ;
  v = vd(d,grp1totaln,grp2totaln) ;
  frm.d.value =  round(d,4);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

function calc28()
{
  var frm = document.frm ;
  fvalue = parseFloat(frm.fvalue.value) ;
  var mean    = new Array() ;
  var n       = new Array() ;
  mean[0] = parseFloat(frm.mean0.value) ;
  mean[1] = parseFloat(frm.mean1.value) ;
  mean[2] = parseFloat(frm.mean2.value) ;
  mean[3] = parseFloat(frm.mean3.value) ;
  mean[4] = parseFloat(frm.mean4.value) ;
  mean[5] = parseFloat(frm.mean5.value) ;
  mean[6] = parseFloat(frm.mean6.value) ;
  mean[7] = parseFloat(frm.mean7.value) ;
  n[0] = parseFloat(frm.n0.value) ;
  n[1] = parseFloat(frm.n1.value) ;
  n[2] = parseFloat(frm.n2.value) ;
  n[3] = parseFloat(frm.n3.value) ;
  n[4] = parseFloat(frm.n4.value) ;
  n[5] = parseFloat(frm.n5.value) ;
  n[6] = parseFloat(frm.n6.value) ;
  n[7] = parseFloat(frm.n7.value) ;

  var rowtype = new Array() ;
  for (i=0;i<=2;i+=1)
  {
    if (frm.rowtype0[i].checked) {rowtype[0]=i ;}
    if (frm.rowtype1[i].checked) {rowtype[1]=i ;}
    if (frm.rowtype2[i].checked) {rowtype[2]=i ;}
    if (frm.rowtype3[i].checked) {rowtype[3]=i ;}
    if (frm.rowtype4[i].checked) {rowtype[4]=i ;}
    if (frm.rowtype5[i].checked) {rowtype[5]=i ;}
    if (frm.rowtype6[i].checked) {rowtype[6]=i ;}
    if (frm.rowtype7[i].checked) {rowtype[7]=i ;}
  }

  totaln = 0 ;
  meanXn = 0 ;
  mean2Xn = 0 ;
  grp1nm = 0 ;
  grp2nm = 0 ;
  grp1n = 0 ;
  grp2n = 0 ;
  k = 0 ;
  k1 = 0 ;
  k2 = 0 ;
  for (i=0;i<=7;i=i+1)
  {
   if (mean[i]!="" && n[i]>0)
   {
   k       += 1 ;
   totaln  += n[i] ;
   meanXn  += mean[i]*n[i] ;
   mean2Xn += mean[i]*mean[i]*n[i] ;
     if (rowtype[i]==0)
     {
     grp1nm += mean[i]*n[i] ;
     grp1n  += n[i] ;
     k1     += 1 ;
     }
     else if (rowtype[i]==1)
     {
     grp2nm += mean[i]*n[i] ;
     grp2n  += n[i] ;
     k2     += 1 ;
     }
   }
  }
  msb =  (mean2Xn - (meanXn*meanXn)/totaln)/ (k-1) ;
  sd_pooled = Math.sqrt(msb/fvalue) ;
  grp1m = grp1nm/grp1n ;
  grp2m = grp2nm/grp2n ;
  d = (grp1m - grp2m)/sd_pooled ;
  v = vd(d,grp1n,grp2n) ;
  frm.d.value =  round(d,4);
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}


/* Compute d based on ANCOVA */
function calc29()
{
  var frm = document.frm ;
  var grp1m = parseFloat(frm.grp1m.value) ;
  var grp1n = parseFloat(frm.grp1n.value) ;
  var grp2m = parseFloat(frm.grp2m.value) ;
  var grp2n = parseFloat(frm.grp2n.value) ;
  var mserror = parseFloat(frm.mserror.value) ;
  var r       = parseFloat(frm.r.value) ;
  totaln = grp1n + grp2n ;
  dferror = totaln - 2 ;
  sd_pooled = Math.sqrt((mserror/(1-r*r))*((dferror-1)/(dferror-2))) ;
  d = (grp1m-grp2m) / sd_pooled ;
  v = vd(d,grp1n,grp2n) ;
  frm.d.value = round(d,4) ;
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}

/* Compute d based on 2-way ANOVA */
function calc30()
{
  var frm = document.frm ;
  var fa      = parseFloat(frm.fa.value) ;
  var fb      = parseFloat(frm.fb.value) ;
  var fab     = parseFloat(frm.fab.value) ;
  var dfb     = parseFloat(frm.dfb.value) ;
  var mserror = parseFloat(frm.mserror.value) ;
  var grp1m = parseFloat(frm.grp1m.value) ;
  var grp1n = parseFloat(frm.grp1n.value) ;
  var grp2m = parseFloat(frm.grp2m.value) ;
  var grp2n = parseFloat(frm.grp2n.value) ;
  dftotal = (grp1n + grp2n) - 1 ;
  dfa = 1 ;
  dfab = dfa*dfb ;
  dferror = dftotal - (dfa+dfb+dfab) ;
  ssa = fa*mserror*dfa ;
  ssb = fb*mserror*dfb ;
  ssab = fab*mserror*dfab ;
  sserror = mserror*dferror ;
  // determine which radio button is checked
  if (frm.sdtype[0].checked)
  {
    sd_pooled = Math.sqrt((ssab+sserror)/(dfab+dferror)) ;
  }
  else if (frm.sdtype[1].checked)
  {
   sd_pooled = Math.sqrt((ssb+ssab+sserror)/(dfb+dfab+dferror)) ;
  }
  d = (grp1m-grp2m) / sd_pooled ;
  v = vd(d,grp1n,grp2n) ;
  frm.d.value = d ;
  frm.v.value = round(v,4);
  frm.lower.value = round(lower_d(d,v),4) ;
  frm.upper.value = round(upper_d(d,v),4) ;
}


/* Compute OR/RR based on 2 by 2 frequency data */
function calc_or1()
{
  var frm = document.frm ;
  var a = parseFloat(frm.a.value) ;
  var b = parseFloat(frm.b.value) ;
  var c = parseFloat(frm.c.value) ;
  var d = parseFloat(frm.d.value) ;
  oddsratio = (a*d)/(b*c) ;
  riskratio = (a/(a+b))/(c/(c+d)) ;
  logged_or  = Math.log(oddsratio) ;
  logged_rr  = Math.log(riskratio) ;
  v_or = 1/a + 1/b + 1/c + 1/d ;
  v_rr = ((b/a)/(a+b)) + ((d/c)/(c+d)) ;
  frm.oddsratio.value =  round(oddsratio,4);
  frm.riskratio.value =  round(riskratio,4);
  frm.logged_or.value =  round(logged_or,4);
  frm.logged_rr.value =  round(logged_rr,4);
  frm.v_or.value =  round(v_or,4);
  frm.v_rr.value =  round(v_rr,4);
  frm.lower_or.value =  round(Math.exp(lower_d(logged_or,v_or)),4);
  frm.upper_or.value =  round(Math.exp(upper_d(logged_or,v_or)),4);
  frm.lower_rr.value =  round(Math.exp(lower_d(logged_rr,v_rr)),4);
  frm.upper_rr.value =  round(Math.exp(upper_d(logged_rr,v_rr)),4);
}


/* Compute OR/RR based on 2 by 2 proportion data */
function calc_or2()
{
  var frm = document.frm ;
  var grp1ps = parseFloat(frm.grp1ps.value) ;
  var grp2ps = parseFloat(frm.grp2ps.value) ;
  var grp1n  = parseFloat(frm.grp1n.value) ;
  var grp2n  = parseFloat(frm.grp2n.value) ;
  a = grp1ps*grp1n ;
  b = (1-grp1ps)*grp1n ;
  c = grp2ps*grp2n ;
  d = (1-grp2ps)*grp2n ;
  oddsratio = (a*d)/(b*c) ;
  riskratio = (a/(a+b))/(c/(c+d)) ;
  logged_or  = Math.log(oddsratio) ;
  logged_rr  = Math.log(riskratio) ;
  v_or = 1/a + 1/b + 1/c + 1/d ;
  v_rr = ((b/a)/(a+b)) + ((d/c)/(c+d)) ;
  frm.oddsratio.value =  round(oddsratio,4);
  frm.riskratio.value =  round(riskratio,4);
  frm.logged_or.value =  round(logged_or,4);
  frm.logged_rr.value =  round(logged_rr,4);
  frm.v_or.value =  round(v_or,4);
  frm.v_rr.value =  round(v_rr,4);
  frm.lower_or.value =  round(Math.exp(lower_d(logged_or,v_or)),4);
  frm.upper_or.value =  round(Math.exp(upper_d(logged_or,v_or)),4);
  frm.lower_rr.value =  round(Math.exp(lower_d(logged_rr,v_rr)),4);
  frm.upper_rr.value =  round(Math.exp(upper_d(logged_rr,v_rr)),4);
}

/* Compute OR/RR based on Phi and Marginal Values */
function calc_or3()
{
  var frm = document.frm ;
  var grp1n  = parseFloat(frm.grp1n.value) ;
  var grp2n  = parseFloat(frm.grp2n.value) ;
  var r      = parseFloat(frm.r.value) ;
  var pevent = parseFloat(frm.pevent.value) ;
  totaln = grp1n + grp2n ;
  ptx = grp1n/totaln ;
  a = totaln *  ((ptx*pevent)+ r*Math.sqrt(ptx*pevent*(1-ptx)*(1-pevent))) ;
  b = grp1n - a ;
  c = totaln*pevent - a ;
  d = totaln - (a+b+c) ;
  oddsratio = (a*d)/(b*c) ;
  riskratio = (a/(a+b))/(c/(c+d)) ;
  logged_or  = Math.log(oddsratio) ;
  logged_rr  = Math.log(riskratio) ;
  v_or = 1/a + 1/b + 1/c + 1/d ;
  v_rr = ((b/a)/(a+b)) + ((d/c)/(c+d)) ;
  frm.oddsratio.value =  round(oddsratio,4);
  frm.riskratio.value =  round(riskratio,4);
  frm.logged_or.value =  round(logged_or,4);
  frm.logged_rr.value =  round(logged_rr,4);
  frm.v_or.value =  round(v_or,4);
  frm.v_rr.value =  round(v_rr,4);
  frm.lower_or.value =  round(Math.exp(lower_d(logged_or,v_or)),4);
  frm.upper_or.value =  round(Math.exp(upper_d(logged_or,v_or)),4);
  frm.lower_rr.value =  round(Math.exp(lower_d(logged_rr,v_rr)),4);
  frm.upper_rr.value =  round(Math.exp(upper_d(logged_rr,v_rr)),4);
}

/* Compute OR/RR based on Chi-square and Marginal Values */
function calc_or4()
{
  var frm = document.frm ;
  var grp1n  = parseFloat(frm.grp1n.value) ;
  var grp2n  = parseFloat(frm.grp2n.value) ;
  var chisq      = parseFloat(frm.chisq.value) ;
  var pevent = parseFloat(frm.pevent.value) ;
  totaln = grp1n + grp2n ;
  ptx = grp1n/totaln ;
  a = totaln*((ptx*pevent)+Math.sqrt((chisq*ptx*pevent*(1-ptx)*(1-pevent))/totaln)) ;
  b = grp1n - a ;
  c = totaln*pevent - a ;
  d = totaln - (a+b+c) ;
  oddsratio = (a*d)/(b*c) ;
  riskratio = (a/(a+b))/(c/(c+d)) ;
  logged_or  = Math.log(oddsratio) ;
  logged_rr  = Math.log(riskratio) ;
  v_or = 1/a + 1/b + 1/c + 1/d ;
  v_rr = ((b/a)/(a+b)) + ((d/c)/(c+d)) ;
  frm.oddsratio.value =  round(oddsratio,4);
  frm.riskratio.value =  round(riskratio,4);
  frm.logged_or.value =  round(logged_or,4);
  frm.logged_rr.value =  round(logged_rr,4);
  frm.v_or.value =  round(v_or,4);
  frm.v_rr.value =  round(v_rr,4);
  frm.lower_or.value =  round(Math.exp(lower_d(logged_or,v_or)),4);
  frm.upper_or.value =  round(Math.exp(upper_d(logged_or,v_or)),4);
  frm.lower_rr.value =  round(Math.exp(lower_d(logged_rr,v_rr)),4);
  frm.upper_rr.value =  round(Math.exp(upper_d(logged_rr,v_rr)),4);
}


/* Compute OR/RR based on Chi-square and Marginal Values */
function calc_or5()
{
  var frm = document.frm ;
  var d  = parseFloat(frm.d.value) ;
  var v  = parseFloat(frm.v.value) ;
  logged_or1 =  Math.PI/Math.sqrt(3) * d ;
  v_or1 = Math.pow(Math.PI,2)/3 * v ;
  logged_or2 = d*1.65 ;
  v_or2 = v/.367 ;
  oddsratio1 = Math.exp(logged_or1) ;
  oddsratio2 = Math.exp(logged_or2) ;
  frm.oddsratio1.value =  round(oddsratio1,4);
  frm.oddsratio2.value =  round(oddsratio2,4);
  frm.logged_or1.value =  round(logged_or1,4);
  frm.logged_or2.value =  round(logged_or2,4);
  frm.v_or1.value =  round(v_or1,4);
  frm.v_or2.value =  round(v_or2,4);
  frm.lower_or1.value =  round(Math.exp(lower_d(logged_or1,v_or1)),4);
  frm.upper_or1.value =  round(Math.exp(upper_d(logged_or1,v_or1)),4);
  frm.lower_or2.value =  round(Math.exp(lower_d(logged_or2,v_or2)),4);
  frm.upper_or2.value =  round(Math.exp(upper_d(logged_or2,v_or2)),4);
}


/* Correlation (r) from k by j frequency table */
function calc_r1()
{
  var frm = document.frm ;
  var freq = new Array(9) ;
  for (i=0; i<9; i++ )
  {
    freq[i] = new Array(9) ;
  }
  //for (k=0; k<9; k++ )
  //{
  //  for (j=0; j<0; j++ )
  //  {
  //  myFormVar = "frm.k" + k + "j" + j  ;
  //  freq[k][j] = parseFloast(myFormVar.value) ;
  //  }
  //}
  freq[0][0]  = parseFloat(frm.k0j0.value) ;  freq[0][1]  = parseFloat(frm.k0j1.value) ;
  freq[0][2]  = parseFloat(frm.k0j2.value) ;  freq[0][3]  = parseFloat(frm.k0j3.value) ;
  freq[0][4]  = parseFloat(frm.k0j4.value) ;  freq[0][5]  = parseFloat(frm.k0j5.value) ;
  freq[0][6]  = parseFloat(frm.k0j6.value) ;  freq[0][6]  = parseFloat(frm.k0j7.value) ;
  freq[0][8]  = parseFloat(frm.k0j8.value) ;  freq[1][0]  = parseFloat(frm.k1j0.value) ;
  freq[1][1]  = parseFloat(frm.k1j1.value) ;  freq[1][2]  = parseFloat(frm.k1j2.value) ;
  freq[1][3]  = parseFloat(frm.k1j3.value) ;  freq[1][4]  = parseFloat(frm.k1j4.value) ;
  freq[1][5]  = parseFloat(frm.k1j5.value) ;  freq[1][6]  = parseFloat(frm.k1j6.value) ;
  freq[1][6]  = parseFloat(frm.k1j7.value) ;  freq[1][8]  = parseFloat(frm.k1j8.value) ;
  freq[2][0]  = parseFloat(frm.k2j0.value) ;  freq[2][1]  = parseFloat(frm.k2j1.value) ;
  freq[2][2]  = parseFloat(frm.k2j2.value) ;  freq[2][3]  = parseFloat(frm.k2j3.value) ;
  freq[2][4]  = parseFloat(frm.k2j4.value) ;  freq[2][5]  = parseFloat(frm.k2j5.value) ;
  freq[2][6]  = parseFloat(frm.k2j6.value) ;  freq[2][6]  = parseFloat(frm.k2j7.value) ;
  freq[2][8]  = parseFloat(frm.k2j8.value) ;  freq[3][0]  = parseFloat(frm.k3j0.value) ;
  freq[3][1]  = parseFloat(frm.k3j1.value) ;  freq[3][2]  = parseFloat(frm.k3j2.value) ;
  freq[3][3]  = parseFloat(frm.k3j3.value) ;  freq[3][4]  = parseFloat(frm.k3j4.value) ;
  freq[3][5]  = parseFloat(frm.k3j5.value) ;  freq[3][6]  = parseFloat(frm.k3j6.value) ;
  freq[3][6]  = parseFloat(frm.k3j7.value) ;  freq[3][8]  = parseFloat(frm.k3j8.value) ;
  freq[4][0]  = parseFloat(frm.k4j0.value) ;  freq[4][1]  = parseFloat(frm.k4j1.value) ;
  freq[4][2]  = parseFloat(frm.k4j2.value) ;  freq[4][3]  = parseFloat(frm.k4j3.value) ;
  freq[4][4]  = parseFloat(frm.k4j4.value) ;  freq[4][5]  = parseFloat(frm.k4j5.value) ;
  freq[4][6]  = parseFloat(frm.k4j6.value) ;  freq[4][6]  = parseFloat(frm.k4j7.value) ;
  freq[4][8]  = parseFloat(frm.k4j8.value) ;  freq[5][0]  = parseFloat(frm.k5j0.value) ;
  freq[5][1]  = parseFloat(frm.k5j1.value) ;  freq[5][2]  = parseFloat(frm.k5j2.value) ;
  freq[5][3]  = parseFloat(frm.k5j3.value) ;  freq[5][4]  = parseFloat(frm.k5j4.value) ;
  freq[5][5]  = parseFloat(frm.k5j5.value) ;  freq[5][6]  = parseFloat(frm.k5j6.value) ;
  freq[5][6]  = parseFloat(frm.k5j7.value) ;  freq[5][8]  = parseFloat(frm.k5j8.value) ;
  freq[6][0]  = parseFloat(frm.k6j0.value) ;  freq[6][1]  = parseFloat(frm.k6j1.value) ;
  freq[6][2]  = parseFloat(frm.k6j2.value) ;  freq[6][3]  = parseFloat(frm.k6j3.value) ;
  freq[6][4]  = parseFloat(frm.k6j4.value) ;  freq[6][5]  = parseFloat(frm.k6j5.value) ;
  freq[6][6]  = parseFloat(frm.k6j6.value) ;  freq[6][6]  = parseFloat(frm.k6j7.value) ;
  freq[6][8]  = parseFloat(frm.k6j8.value) ;  freq[7][0]  = parseFloat(frm.k7j0.value) ;
  freq[7][1]  = parseFloat(frm.k7j1.value) ;  freq[7][2]  = parseFloat(frm.k7j2.value) ;
  freq[7][3]  = parseFloat(frm.k7j3.value) ;  freq[7][4]  = parseFloat(frm.k7j4.value) ;
  freq[7][5]  = parseFloat(frm.k7j5.value) ;  freq[7][6]  = parseFloat(frm.k7j6.value) ;
  freq[7][6]  = parseFloat(frm.k7j7.value) ;  freq[7][8]  = parseFloat(frm.k7j8.value) ;
  freq[8][0]  = parseFloat(frm.k8j0.value) ;  freq[8][1]  = parseFloat(frm.k8j1.value) ;
  freq[8][2]  = parseFloat(frm.k8j2.value) ;  freq[8][3]  = parseFloat(frm.k8j3.value) ;
  freq[8][4]  = parseFloat(frm.k8j4.value) ;  freq[8][5]  = parseFloat(frm.k8j5.value) ;
  freq[8][6]  = parseFloat(frm.k8j6.value) ;  freq[8][6]  = parseFloat(frm.k8j7.value) ;
  freq[8][8]  = parseFloat(frm.k8j8.value) ;

  fkj = 0 ; fk = 0 ; fj = 0 ; fk2 = 0 ; fj2 = 0 ; n = 0 ;
  for (k=0; k<=8; k++)
  {
   for (j=0; j<=8; j++)
   {
     if (freq[k][j]>0)
     {
     fkj += freq[k][j] * k * j ;
     fk  += freq[k][j] * k ;
     fj  += freq[k][j] * j ;
     fk2 += freq[k][j] * k*k ;
     fj2 += freq[k][j] * j*j ;
     n   += freq[k][j] ;
     }
   }
  }
  r = (n*fkj-fk*fj)/Math.sqrt((n*fk2 - fk*fk)*(n*fj2-fj*fj))
  v = 1/(n-3) ;
  zr = Zr(r) ;
  frm.r.value = round(r,4);
  frm.lowerr.value = round(InvZr(lower_d(zr,v)),4) ;
  frm.upperr.value = round(InvZr(upper_d(zr,v)),4) ;
  frm.v.value = round(v,4);
  frm.zr.value = round(zr,4) ;
  frm.lowerz.value = round(lower_d(zr,v),4) ;
  frm.upperz.value = round(upper_d(zr,v),4) ;
}

/* Correlation (r) from r */
function calc_r2()
{
  var frm = document.frm ;
  var r   = parseFloat(frm.rin.value) ;
  var n   = parseFloat(frm.n.value) ;
  v = 1/(n-3) ;
  zr = Zr(r) ;
  frm.rout.value = round(r,4);
  frm.lowerr.value = round(InvZr(lower_d(zr,v)),4) ;
  frm.upperr.value = round(InvZr(upper_d(zr,v)),4) ;
  frm.v.value = round(v,4);
  frm.zr.value = round(zr,4) ;
  frm.lowerz.value = round(lower_d(zr,v),4) ;
  frm.upperz.value = round(upper_d(zr,v),4) ;
}

/* Correlation (r) from Means and SDs */
function calc_r3()
{
  var frm = document.frm ;
  var grp1m = parseFloat(frm.grp1m.value) ;
  var grp1sd = parseFloat(frm.grp1sd.value) ;
  var grp1n = parseFloat(frm.grp1n.value) ;
  var grp2m = parseFloat(frm.grp2m.value) ;
  var grp2sd = parseFloat(frm.grp2sd.value) ;
  var grp2n = parseFloat(frm.grp2n.value) ;
  sd_pooled = Math.sqrt((grp1sd*grp1sd*(grp1n-1)+grp2sd*grp2sd*(grp2n-1))/(grp1n+grp2n-2)) ;
  d = (grp1m-grp2m) / sd_pooled ;
  vd = vd(d,grp1n,grp2n);
  p = grp1n/(grp1n+grp2n) ;
  r = d / Math.sqrt(d*d + 1/(p*(1-p))) ;
  v = vd/(vd + 1/(p*(1-p))) ;
  zr = Zr(r) ;
  frm.r.value = round(r,4);
  frm.lowerr.value = round(InvZr(lower_d(zr,v)),4) ;
  frm.upperr.value = round(InvZr(upper_d(zr,v)),4) ;
  frm.v.value = round(v,4);
  frm.zr.value = round(zr,4) ;
  frm.lowerz.value = round(lower_d(zr,v),4) ;
  frm.upperz.value = round(upper_d(zr,v),4) ;
}

/* Correlation (r) from 2 by 2 */
function calc_r4()
{
  var frm = document.frm ;
  var a = parseFloat(frm.a.value) ;
  var b = parseFloat(frm.b.value) ;
  var c = parseFloat(frm.c.value) ;
  var d = parseFloat(frm.d.value) ;
  r = (a*d - b*c) / Math.sqrt((a+b)*(c+d)*(a+c)*(b+d)) ;
  n = a+b+c+d ;
  // base v on a rescaled v from a logged odds-ratio (this is more precise than n-3
  logged_or  = Math.log((a*d)/(b*c)) ;
  v_or = 1/a + 1/b + 1/c + 1/d ;
  zr = Zr(r) ;
  v = ((zr*zr)/(logged_or*logged_or)) * v_or
  frm.r.value = round(r,4);
  frm.lowerr.value = round(InvZr(lower_d(zr,v)),4) ;
  frm.upperr.value = round(InvZr(upper_d(zr,v)),4) ;
  frm.v.value = round(v,4);
  frm.zr.value = round(zr,4) ;
  frm.lowerz.value = round(lower_d(zr,v),4) ;
  frm.upperz.value = round(upper_d(zr,v),4) ;
}

/* chi-square and N */
function calc_r5()
{
  var frm = document.frm ;
  var chisq   = parseFloat(frm.chisq.value) ;
  var n   = parseFloat(frm.n.value) ;
  r = Math.sqrt(chisq/n) ;
  zr = Zr(r) ;
  v = (zr*zr)/chisq ;
  frm.r.value = round(r,4) ;
  frm.lowerr.value = round(InvZr(lower_d(zr,v)),4) ;
  frm.upperr.value = round(InvZr(upper_d(zr,v)),4) ;
  frm.v.value = round(v,4);
  frm.zr.value = round(zr,4) ;
  frm.lowerz.value = round(lower_d(zr,v),4) ;
  frm.upperz.value = round(upper_d(zr,v),4) ;
}

/* t-test and N */
function calc_r6()
{
  var frm = document.frm ;
  var t   = parseFloat(frm.t.value) ;
  var n   = parseFloat(frm.n.value) ;
  r = t / Math.sqrt( (t*t) + n - 2) ;
  zr = Zr(r) ;
  v = 1/(n-3) ;
  frm.r.value = round(r,4);
  frm.lowerr.value = round(InvZr(lower_d(zr,v)),4) ;
  frm.upperr.value = round(InvZr(upper_d(zr,v)),4) ;
  frm.v.value = round(v,4);
  frm.zr.value = round(zr,4) ;
  frm.lowerz.value = round(lower_d(zr,v),4) ;
  frm.upperz.value = round(upper_d(zr,v),4) ;
}

/* t-test p-value and N */
function calc_r7()
{
  var frm = document.frm ;
  var poft = parseFloat(frm.poft.value) ;
  var n = parseFloat(frm.n.value) ;
  t = tfromp(n-2,poft) ;
  r = t / Math.sqrt( (t*t) + n - 2) ;
  zr = Zr(r) ;
  v = 1/(n-3) ;
  frm.r.value = round(r,4);
  frm.lowerr.value = round(InvZr(lower_d(zr,v)),4) ;
  frm.upperr.value = round(InvZr(upper_d(zr,v)),4) ;
  frm.v.value = round(v,4);
  frm.zr.value = round(zr,4) ;
  frm.lowerz.value = round(lower_d(zr,v),4) ;
  frm.upperz.value = round(upper_d(zr,v),4) ;
}


/* Odeh and Evans (1974) Inverse Normal Distribution Function */
/* See Brophy 1985 */
function zfromp(pvalue)
{
  var y = Math.sqrt( -2 * Math.log(pvalue/2)) ;
  var z = y - ((((0.0000453642210148*y + 0.0204231210245)*y +
        0.342242088547)*y+1)*y + 0.322232431088)/
        ((((0.0038560700634*y+0.10353775285)*y + 0.531103462366)*y
        + 0.588581570495)*y+0.099348462606) ;
  return z ;
}


/* Hill 1970 Students t-quantiles; Algorithm 396 */
function tfromp(df,pvalue)
{
    var halfpi = Math.PI/2 ;
    var a = 1/(df-.5) ;
    var b = 48/(a*a) ;
    var c = ((20700*a/b-98)*a-16)*a+96.36 ;
    var d = ((94.5/(b+c)-3)/b+1)*Math.sqrt(a*halfpi)*df ;
    var x = d * pvalue ;
    var y = Math.pow(x,2/df) ;
    if (y>(.05+a)) {
	var x = zfromp(pvalue)*-1 ;
        var y = Math.pow(x,2) ;
        var c = (((.05*d*x-5)*x-7)*x-2)*x+b+c ;
        var y = (((((0.4*y+6.3)*y+36)*y+94.5)/c-y-3)/b+1)*x ;
        var y = a * Math.pow(y,2) ;
        if (y>.002) { t = Math.sqrt((Math.exp(y)-1)*df)  } else  { t = Math.sqrt((.5 * Math.pow(y,2) + y)*df)  }
    }
    else  {
        var y = ((1/(((df+6)/(df*y)-0.089*d-0.822)*(df+2)*3) +
                  0.5/(df+4))*y-1)*(df+1)/(df+2)+1/y ;
        var t = Math.sqrt(df*y)
    }
	return t ;
}

