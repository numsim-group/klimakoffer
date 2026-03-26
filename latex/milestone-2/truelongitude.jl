
"""
calc_lambda()
Solve the orbital equation for lambda as a function of time with
a fourth-order Runge-Kutta method  
"""
function calc_lambda(dt, ecc, per, nt)

    eccfac = 1.0 - ecc^2
    rzero  = (2.0*pi)/eccfac^1.5
  
    lambda = zeros(Float64,nt+1)
  
    for n in 2:nt+1
      nu = lambda[n-1] - per
      t1 = dt*(rzero*(1.0 - ecc * cos(nu))^2)
      t2 = dt*(rzero*(1.0 - ecc * cos(nu+0.5*t1))^2)
      t3 = dt*(rzero*(1.0 - ecc * cos(nu+0.5*t2))^2)
      t4 = dt*(rzero*(1.0 - ecc * cos(nu + t3))^2)
      lambda[n] = lambda[n-1] + (t1 + 2.0*t2 + 2.0*t3 + t4)/6.0
    end
  
    return lambda
  end

  dt = 1.0 / 48
  ecc = 0.016740  
  per = 1.783037 
  nt = 48

  lambda = calc_lambda(dt, ecc, per, nt)

display(lambda)