# upperEV.R upper extreme value copula

#' Joint cdf of the upper extreme value copula of the GGEE copula
#'
#' Joint cdf of the upper extreme value copula of the GGEE copula
#' @param u, v: inputs for the extreme value copula
#' @param b: the parameter of the extreme value copula
#' @keywords upper extreme value copula
#' @export
#' @examples
#' pUEV_GGEE_COP(0.3, 0.4, 1.2)
pUEV_GGEE_COP <- function(u, v, b)
{
  if (u == v)
  {
    out <- u^(1 + b)
  } else
  {
    w1 <- -log(u)
    w2 <- -log(v)
    out <- exp(-(w1^(1 + 1/b) - w2^(1 + 1/b))/(w1^(1/b) - w2^(1/b)))
  }
  return(out)
}

pow <- function(x, y){ x^y }

#' Density function of the upper extreme value copula for the GGEE copula
#'
#' Density function of the upper extreme value copula for the GGEE copula
#' @param u, v: inputs for the extreme value copula
#' @param b: the parameter of the extreme value copula
#' @keywords upper extreme value copula
#' @export
#' @examples
#' dUEV_GGEE_COP(0.3, 0.4, 1.2)
dUEV_GGEE_COP <- function(u, v, b)
{
  if (u == v)
  {
    out <- (1/12) * (3 * v^b * b^3 * log(v) + 6 * v^b * b^2 * log(v) + 2 * 
      v^b * b^2 + 3 * b * log(v) * v^b - 2 * v^b)/(v * log(v) * b)
  } else
  {
    t1 <- log(u)
    t2 <- 1/b
    t3 <- pow(-t1, t2)
    t4 <- log(v)
    t5 <- pow(-t4, t2)
    t6 <- t3 - t5
    t7 <- 1/t6
    t9 <- pow(u, t7 * t5)
    t11 <- pow(v, -t7 * t3)
    t13 <- 2 * t2
    t14 <- pow(-t1, t13)
    t16 <- (-2 + b) * t2
    t17 <- pow(-t4, -t16)
    t21 <- 3 * t2
    t22 <- pow(-t4, t21)
    t23 <- t3 * t22
    t25 <- pow(-t1, t21)
    t26 <- t25 * t5
    t28 <- pow(-t4, t13)
    t29 <- t14 * t28
    t32 <- (-3 + b) * t2
    t33 <- pow(-t4, -t32)
    t34 <- t3 * t33
    t36 <- pow(-t1, -t32)
    t37 <- t5 * t36
    t39 <- pow(-t1, -t16)
    t44 <- (-1 + b) * t2
    t45 <- pow(-t4, -t44)
    t46 <- t25 * t45
    t47 <- pow(-t1, -t44)
    t48 <- t22 * t47
    t49 <- -2 * t14 * t17 * b + t23 * b + t26 * b - 2 * t29 + t34 * b + t37 * 
      b - 2 * t28 * t39 * b - t46 - t48 + t34 + t37
    t51 <- (2 + b) * t2
    t52 <- pow(-t1, t51)
    t53 <- t52 * t17
    t54 <- pow(-t4, t51)
    t55 <- t39 * t54
    t58 <- (3 + b) * t2
    t59 <- pow(-t4, t58)
    t66 <- b * b
    t68 <- pow(-t1, t58)
    t75 <- t53 + t55 + t55 * b - t47 * t59 * b - 2 * t29 * b + t46 * b + t48 * 
      b + t26 * t66 - t68 * t45 * b + t53 * b - 2 * t29 * t66 + t23 * t66
    t79 <- t6 * t6
    t80 <- t79 * t79
    out <- -t9 * t11 * (t49 + t75)/t66/t80
  }
  return(out)
}

A1_UEV_PPPP_COP <- function(be,a,b)
{
  be*(a/(a+b)/(be-b) - a*b/((a+b)^2)/(a-b+be) + ((a/(a+b))^2)/(2*b-be))
}

A2_UEV_PPPP_COP <- function(be,a,b)
{
  be*(a*b/((a+b)^2)/(a-b+be) - b/(a+b)/(a+be) + ((a/(a+b))^2)/(2*a+be))
}

A2_UEV_PPPP_COP <- function(be,a,b)
{
  be*(1/be - b/(a+b)/(a+be) - a/(a+b)/(be-b))
}

A_UEV_PPPP_COP <- function(w1,w2,be,a,b)
{
  wm <- min(w1, w2)
  wp <- max(w1, w2)
  wm+wp-(a+be)*(b-be)/(a*b)*(A1_UEV_PPPP_COP(be,a,b)*(wm^(b/be))*(wp^(1-b/be)) + A2_UEV_PPPP_COP(be,a,b)*(wm^(1+a/be))*(wp^(-a/be)) + A3_UEV_PPPP_COP(be,a,b)*wm )
}

#' Joint cdf of the upper extreme value copula of the PPPP copula
#'
#' Joint cdf of the upper extreme value copula of the PPPP copula
#' @param u, v: inputs for the extreme value copula
#' @param be,a,b: the parameters of the extreme value copula
#' @keywords upper extreme value copula
#' @export
#' @examples
#' pUEV_PPPP_COP(0.2, 0.4, 1.2,1,1)
pUEV_PPPP_COP <- function(u, v, be, a, b)
{
  if (b <= be)
  {
    out <- u*v
  } else
  {
    w1 <- -log(u)
    w2 <- -log(v)
    out <- exp(-A_UEV_PPPP_COP(w1,w2,be,a,b))
  }
  return(out)
}


#' Density function of the upper extreme value copula for the PPPP copula
#'
#' Density function of the upper extreme value copula for the PPPP copula
#' @param u, v: inputs for the extreme value copula
#' @param be,a,b: the parameter of the extreme value copula
#' @keywords upper extreme value copula
#' @export
#' @examples
#' dUEV_PPPP_COP(0.3, 0.4, 0.6, 1, 1)
dUEV_PPPP_COP <- function(u1,u2,be,a,b)
{
  if(u1!=u2)
  {
    if(u2<u1)
    {
      u0 <- u1
      u1 <- u2
      u2 <- u0
    }
    
    t1 = a+be;
    t2 = b-be;
    t3 = t1*t2;
    t4 = 1/a;
    t6 = 1/b;
    t7 = a+b;
    t8 = 1/t7;
    t11 = -a*t8/t2;
    t13 = t7*t7;
    t14 = 1/t13;
    t18 = a*b*t14/(a-b+be);
    t19 = a*a;
    t25 = t11-t18+t19*t14/(2.0*b-be);
    t26 = log(u2);
    t27 = 1/be;
    t28 = b*t27;
    t29 = pow(-t26,1.0*t28);
    t30 = t25*t29;
    t31 = 1/u2;
    t34 = 1/t26;
    t35 = log(u1);
    t36 = 1.0-t28;
    t37 = pow(-t35,1.0*t36);
    t39 = 1/u1;
    t41 = 1/t35;
    t47 = b*t8/t1;
    t48 = b*b;
    t54 = t18-t47+t48*t14/(2.0*a+be);
    t55 = a*t27;
    t56 = 1.0+t55;
    t57 = pow(-t26,1.0*t56);
    t58 = t54*t57;
    t59 = t56*t31;
    t61 = pow(-t35,-1.0*t55);
    t62 = t34*t61;
    t64 = a*t39*t41;
    t69 = t4*t6;
    t70 = be*t25;
    t73 = be*t54;
    t77 = be*(t27-t47-t11);
    t83 = exp(t26+t35+t3*t69*(t70*t29*t37+t73*t57*t61-t77*t26));
    out = t3*t4*t6*(t30*b*t31*t34*t37*t36*t39*t41-t58*t59*t62*t64)*t83+(t39+t3*t69*(t70*t29*t37*t36*t39*t41-t58*t61*t64))*(t31+t3*t69*(t30*b*t31*t34*t37+t73*t57*t59*t62-t77*t31))*t83;   
  }else
  {
    t1 = a+be;
    t2 = b-be;
    t3 = t1*t2;
    t4 = 1/a;
    t6 = 1/b;
    t7 = a+b;
    t8 = 1/t7;
    t11 = -a*t8/t2;
    t13 = t7*t7;
    t14 = 1/t13;
    t18 = a*b*t14/(a-b+be);
    t19 = a*a;
    t25 = t11-t18+t19*t14/(2.0*b-be);
    t26 = log(u1);
    t27 = 1/be;
    t28 = b*t27;
    t29 = pow(-t26,1.0*t28);
    t31 = t25*t29*b;
    t32 = u1*u1;
    t34 = t26*t26;
    t36 = 1/t32/t34;
    t37 = 1.0-t28;
    t38 = pow(-t26,1.0*t37);
    t39 = t38*t37;
    t44 = b*t8/t1;
    t45 = b*b;
    t51 = t18-t44+t45*t14/(2.0*a+be);
    t52 = a*t27;
    t53 = 1.0+t52;
    t54 = pow(-t26,1.0*t53);
    t55 = t51*t54;
    t57 = pow(-t26,-1.0*t52);
    t64 = t4*t6;
    t65 = be*t25;
    t68 = be*t51;
    t72 = be*(t27-t44-t11);
    t78 = exp(2.0*t26+t3*t64*(t65*t29*t38+t68*t54*t57-t72*t26));
    t81 = 1/u1;
    t83 = 1/t26;
    t84 = t81*t83;
    out = t3*t4*t6*(t31*t36*t39-t55*t53*t36*t57*a)*t78+(t81+t3*t64*(t65*t29*t39*t84-t55*t57*a*t81*t83))*(t81+t3*t64*(t31*t84*t38+t68*t54*t53*t81*t83*t57-t72*t81))*t78;
  }
  return(out)
}



