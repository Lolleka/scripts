function [ output_args ] = BPPDistribution( x, varargin )
    tau0 = varargin{1};
    DPos = varargin{2};
    DWidth  = varargin{3};
    DType = varargin{4}; % 0 = gauss, 1 = lorentz
    Ampl = varargin{5};
    OmegaL = varargin{6};
    
    
    Distr = 0;
    switch DType
        case 0 % gauss
            lowlim = max([DPos-5*DWidth 0]);
            highlim = DPos+5*DWidth;
            step = DWidth/100;
            E = lowlim:step:highlim;
            Distr = (1/(DWidth*sqrt(2*pi))).*exp(-(E-DPos).*(E-DPos)./(2*DWidth*DWidth));
        case 1 % lorentz
            Distr = 0;
    end
    output_args = x*Distr(1);
    
    %figure(2);
    %plot (E,Distr,'-');
    
    y = zeros(size(x));
    for i = 1:length(x)
        tau = tau0*exp(E./x(i));
        BPP = Distr.*tau./(1+OmegaL*OmegaL*tau.*tau);
        y(i) = trapz(E,BPP);
    end
    
    output_args = Ampl*y;
    
    %
    %BPP(T) = int_0^inf(Distr(E)*BPP(T,E))
end