
function xdot = deriv(x,u, K_si_x,limit)


  K_si_si = .0667;
  K_v = 1;

  si_dot = K_si_x * (u(1) - x(1)) + K_si_si * (u(3) - x(4));
  if limit
      if si_dot > 3
          si_dot = 3;
      elseif si_dot < -3
          si_dot = -3;
      else
          si_dot;
      end
  end
  xdot = [1.68781 * x(3) * sind(x(4));
          1.68781 * x(3) * cosd(x(4));
          K_v * (u(2) - x(3));
          si_dot];
end
