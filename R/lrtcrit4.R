#' Find test statistic value from a given data and corresponding critical value
#'
#' More detailed description
#'
#' @param X1 a real vector
#' @param X2 a real vector
#' @param X3 a real vector
#' @param X4 a real vector
#' @param alpha real number between 0 and 1 called significance level
#'
#' @return a numeric vector
#'
#' @examples
#' X1=c(4.28,2.56,0.83,1.41,2.18,1.47,2.66,-0.73,3.58,1.55,4.2,2.98,1.12,0.22,4.3)
#' X2=c(2.12,-1.03,3.37,6.95,2.42,4.22,5.97,0.84,4.57,-0.16,2.18,2.56)
#' X3=c(4.76,3.85,3.54,3.55,2.73,2.61,3.81,3.61,3.28,3.01,3.67,1.56,4.6,4.74)
#' X4=c(4.7,3.74,5.82,2.03,4.9,4.69,1.45,4.9,2.5,3.33)
#' lrtcrit4(X1,X2,X3,X4,0.05)
#'
#' @export
lrtcrit4<-function(X1,X2,X3,X4,alpha)
{
  fun1<-function(n1,n2,n3,n4,v1,v2,v3,v4)
  {
    X1=rnorm(n1,0,sqrt(v1))
    X2=rnorm(n2,0,sqrt(v2))
    X3=rnorm(n3,0,sqrt(v3))
    X4=rnorm(n4,0,sqrt(v4))
    m1<-mean(X1);m2<-mean(X2);m3<-mean(X3);m4<-mean(X4)
    x<-c(m1, m2, m3, m4)
    s1<-sum((X1-m1)^2)/n1
    s2<-sum((X2-m2)^2)/n2
    s3<-sum((X3-m3)^2)/n3
    s4<-sum((X4-m4)^2)/n4
    s<-c(s1, s2, s3, s4)
    x0=x
    s0<-s
    x8<-x
    repeat
    {
      w1=c(n1/s0[1],n2/s0[2],n3/s0[3],n4/s0[4])
      x6=pava(x0, w1)
      s11<-sum((X1-x6[1])^2)/n1
      s12<-sum((X2-x6[2])^2)/n2
      s13<-sum((X3-x6[3])^2)/n3
      s14<-sum((X4-x6[4])^2)/n4
      v=c(s11, s12, s13, s14)
      s0<-v
      x10=abs(x8-x6)
      p=max(x10[1], x10[2], x10[3], x10[4], na.rm = FALSE)
      if(p<=0.0000001)
      {
        break
      }
      x8<-x6
    }
    p1<-sum(X1^2)/n1;p2<-sum(X2^2)/n2;p3<-sum(X3^2)/n3;p4<-sum(X4^2)/n4
    f1<-p2*p3
    f2<--2*(m3*p2+m2*p3)
    f3<-p2+4*m2*m3+p3
    f4<--2*(m2+m3)
    g1<-p1*p3
    g2<--2*(m3*p1+m1*p3)
    g3<-p3+4*m1*m3+p1
    g4<--2*(m1+m3)
    h1<-p1*p2
    h2<--2*(m1*p2+m2*p1)
    h3<-p1+p2+4*m1*m2
    h4<--2*(m1+m2)
    a1<-f1*p4
    a2<-f2*p4-2*f1*m4
    a3<-f1-2*f2*m4+f3*p4
    a4<-f2-2*f3*m4+f4*p4
    a5<-f3-2*f4*m4+p4
    a6<-f4-2*m4
    b1<-g1*p4
    b2<-g2*p4-2*g1*m4
    b3<-g1-2*g2*m4+g3*p4
    b4<-g2-2*g3*m4+g4*p4
    b5<-g3-2*g4*m4+p4
    b6<-g4-2*m4
    c1<-h1*p4
    c2<-h2*p4-2*m4*h1
    c3<-h1-2*h2*m4+h3*p4
    c4<-h2-2*h3*m4+h4*p4
    c5<-h3-2*h4*m4+p4
    c6<-h4-2*m4
    d1<-h1*p3
    d2<-h2*p3-2*m3*h1
    d3<-h1-2*h2*m3+h3*p3
    d4<-h2-2*h3*m3+h4*p3
    d5<-h3-2*h4*m3+p3
    d6<-h4-2*m3
    #coefficients of the polynomial
    z1<-n1*m1*a1+n2*m2*b1+n3*m3*c1+n4*m4*d1
    z2<-n1*(m1*a2-a1)+n2*(m2*b2-b1)+n3*(m3*c2-c1)+n4*(m4*d2-d1)
    z3<-n1*(m1*a3-a2)+n2*(m2*b3-b2)+n3*(m3*c3-c2)+n4*(m4*d3-d2)
    z4<-n1*(m1*a4-a3)+n2*(m2*b4-b3)+n3*(m3*c4-c3)+n4*(m4*d4-d3)
    z5<-n1*(m1*a5-a4)+n2*(m2*b5-b4)+n3*(m3*c5-c4)+n4*(m4*d5-d4)
    z6<-n1*(m1*a6-a5)+n2*(m2*b6-b5)+n3*(m3*c6-c5)+n4*(m4*d6-d5)
    z7<-n1*(m1-a6)+n2*(m2-b6)+n3*(m3-c6)+n4*(m4-d6)
    z8<--(n1+n2+n3+n4)
    root<-polyroot(c(z1,z2,z3,z4,z5,z6,z7,z8))
    Rert=Re(root)[abs(Im(root))<=0.000001]
    n=length(Rert)
    m7=NULL;m8=NULL;m9=NULL;m10=NULL;m11=NULL
    for(i in 1:n)
    {

      m7[i]=(sum((X1-Rert[i])^2))/n1
      m8[i]=(sum((X2-Rert[i])^2))/n2
      m9[i]=(sum((X3-Rert[i])^2))/n3
      m10[i]=(sum((X4-Rert[i])^2))/n4
      m11[i]=(m7[i]/s0[1])^(n1/2)*(m8[i]/s0[2])^(n2/2)*(m9[i]/s0[3])^(n3/2)*(m10[i]/s0[4])^(n4/2)

    }
    m111=min(m11)
    #neu<-(s0[1])^(n1/2)*(s0[2])^(n2/2)*(s0[3])^(n3/2)*(s0[4])^(n4/2)
    #dno<-m111
    #value=neu/dno
    value=1/m111
    return(value)
  }
  fun2<-function(n1,n2,n3,n4,alpha,v1,v2,v3,v4)
  {
    x<-replicate(5000,fun1(n1,n2,n3,n4,v1,v2,v3,v4))
    y<-sort(x,decreasing=FALSE)
    m=5000*alpha
    c<-y[m]
    return(c)
  }
  fun3<-function(n1,n2,n3,n4,alpha,v1,v2,v3,v4)
  {
    z<-replicate(10,fun2(n1,n2,n3,n4,alpha,v1,v2,v3,v4))
    cr=mean(z)
    return(cr)
  }
  fun4<-function(X1,X2,X3,X4,alpha)
  {
    m1<-mean(X1);m2<-mean(X2);m3<-mean(X3);m4<-mean(X4)
    n1=length(X1);n2=length(X2);n3=length(X3);n4=length(X4)
    x<-c(m1, m2, m3, m4)
    s1<-sum((X1-m1)^2)/n1
    s2<-sum((X2-m2)^2)/n2
    s3<-sum((X3-m3)^2)/n3
    s4<-sum((X4-m4)^2)/n4
    s<-c(s1, s2, s3, s4)
    x0=x
    s0<-s
    x8<-x
    repeat
    {
      w1=c(n1/s0[1],n2/s0[2],n3/s0[3],n4/s0[4])
      x6=pava(x0, w1)
      s11<-sum((X1-x6[1])^2)/n1
      s12<-sum((X2-x6[2])^2)/n2
      s13<-sum((X3-x6[3])^2)/n3
      s14<-sum((X4-x6[4])^2)/n4
      v=c(s11, s12, s13, s14)
      s0<-v
      x10=abs(x8-x6)
      p=max(x10[1], x10[2], x10[3], x10[4], na.rm = FALSE)
      if(p<=0.0000001)
      {
        break
      }
      x8<-x6
    }

    p1<-sum(X1^2)/n1;p2<-sum(X2^2)/n2;p3<-sum(X3^2)/n3;p4<-sum(X4^2)/n4
    f1<-p2*p3
    f2<--2*(m3*p2+m2*p3)
    f3<-p2+4*m2*m3+p3
    f4<--2*(m2+m3)
    g1<-p1*p3
    g2<--2*(m3*p1+m1*p3)
    g3<-p3+4*m1*m3+p1
    g4<--2*(m1+m3)
    h1<-p1*p2
    h2<--2*(m1*p2+m2*p1)
    h3<-p1+p2+4*m1*m2
    h4<--2*(m1+m2)
    a1<-f1*p4
    a2<-f2*p4-2*f1*m4
    a3<-f1-2*f2*m4+f3*p4
    a4<-f2-2*f3*m4+f4*p4
    a5<-f3-2*f4*m4+p4
    a6<-f4-2*m4
    b1<-g1*p4
    b2<-g2*p4-2*g1*m4
    b3<-g1-2*g2*m4+g3*p4
    b4<-g2-2*g3*m4+g4*p4
    b5<-g3-2*g4*m4+p4
    b6<-g4-2*m4
    c1<-h1*p4
    c2<-h2*p4-2*m4*h1
    c3<-h1-2*h2*m4+h3*p4
    c4<-h2-2*h3*m4+h4*p4
    c5<-h3-2*h4*m4+p4
    c6<-h4-2*m4
    d1<-h1*p3
    d2<-h2*p3-2*m3*h1
    d3<-h1-2*h2*m3+h3*p3
    d4<-h2-2*h3*m3+h4*p3
    d5<-h3-2*h4*m3+p3
    d6<-h4-2*m3
    #coefficients of the polynomial
    z1<-n1*m1*a1+n2*m2*b1+n3*m3*c1+n4*m4*d1
    z2<-n1*(m1*a2-a1)+n2*(m2*b2-b1)+n3*(m3*c2-c1)+n4*(m4*d2-d1)
    z3<-n1*(m1*a3-a2)+n2*(m2*b3-b2)+n3*(m3*c3-c2)+n4*(m4*d3-d2)
    z4<-n1*(m1*a4-a3)+n2*(m2*b4-b3)+n3*(m3*c4-c3)+n4*(m4*d4-d3)
    z5<-n1*(m1*a5-a4)+n2*(m2*b5-b4)+n3*(m3*c5-c4)+n4*(m4*d5-d4)
    z6<-n1*(m1*a6-a5)+n2*(m2*b6-b5)+n3*(m3*c6-c5)+n4*(m4*d6-d5)
    z7<-n1*(m1-a6)+n2*(m2-b6)+n3*(m3-c6)+n4*(m4-d6)
    z8<--(n1+n2+n3+n4)
    root<-polyroot(c(z1,z2,z3,z4,z5,z6,z7,z8))
    Rert=Re(root)[abs(Im(root))<=0.01]
    n=length(Rert)
    m7=NULL;m8=NULL;m9=NULL;m10=NULL;m11=NULL
    for(i in 1:n)
    {

      m7[i]=(sum((X1-Rert[i])^2))/n1
      m8[i]=(sum((X2-Rert[i])^2))/n2
      m9[i]=(sum((X3-Rert[i])^2))/n3
      m10[i]=(sum((X4-Rert[i])^2))/n4
      m11[i]=(m7[i]/s0[1])^(n1/2)*(m8[i]/s0[2])^(n2/2)*(m9[i]/s0[3])^(n3/2)*(m10[i]/s0[4])^(n4/2)
    }
    m111=min(m11)
    #neu<-(s0[1])^(n1/2)*(s0[2])^(n2/2)*(s0[3])^(n3/2)*(s0[4])^(n4/2)
    #dno<-m111
    #value=neu/dno
    value=1/m111
    return(value)
  }
  set.seed(39)
  n1=length(X1);n2=length(X2);n3=length(X3);n4=length(X4)
  v1=var(X1);v2=var(X2);v3=var(X3);v4=var(X4)
  statistic_value<-fun4(X1,X2,X3,X4,alpha)
  critical_value<-fun3(n1,n2,n3,n4,alpha,v1,v2,v3,v4)
  result<-c(statistic_value,critical_value )
  return(result)
}
