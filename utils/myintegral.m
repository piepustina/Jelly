function int = myintegral(varargin)
%MYINTEGRAL Implements a numerical integrator. This function is called with
%the same arguments of the integral function of MATLAB.
%int = integral(varargin{:});
% TODO: This function has to be generalized so as to use as many Gaussian
% points as we want
% x_gauss_5 = [-0.906;
% -0.538;
%      0;
%  0.538;
%  0.906];
% w_gauss_5 = [0.2369;
% 0.4786;
% 0.5689;
% 0.4786;
% 0.2369];

x_gauss_10 = [-0.9739;
   -0.8651;
   -0.6794;
   -0.4334;
   -0.1489;
    0.1489;
    0.4334;
    0.6794;
    0.8651;
    0.9739];
w_gauss_10 = [0.0667;
    0.1495;
    0.2191;
    0.2693;
    0.2955;
    0.2955;
    0.2693;
    0.2191;
    0.1495;
    0.0667];

int = int_approx_numerical(varargin{1}, [varargin{2}, varargin{3}], x_gauss_10, w_gauss_10);
end

