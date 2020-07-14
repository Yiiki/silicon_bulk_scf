% LDA
function y=rmc(x)
y=x.^0-x.^0;
a=0.0621814;
b=3.72744;
c=12.9352;
x0=-0.10498;
    for i_x=1:length(x)
        xx=x(i_x);
        if abs(xx)<1e-180
            y(i_x)=0;
        else
            y(i_x)=-a./6.*(c.*((3./(4*pi.*xx)).^(1/6)-x0)-b.*((3./(4*pi.*xx)).^(1/6))*x0)...
                ./(((3./(4*pi.*xx)).^(1/6)-x0).*(((3./(4*pi.*xx)).^(1/3))+b.*((3./(4*pi.*xx)).^(1/6))+c));
        end
    end
end  