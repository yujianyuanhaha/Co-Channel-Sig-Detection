function out = wrap2pi(in)
    out = rem(in,2*pi);
%     out = in;
%     while(jnkisempty(out(out >= pi)))
        out(out>=pi) = out(out >= pi) - 2*pi;
%     end
%     while(jnkisempty(out(out < -pi)))
        out(out < -pi) = out(out < -pi) + 2*pi;
%     end
%     for i=1:length(in)
%         while(out(i) > pi)
%             out(i) = out(i) - 2*pi;
%         end
%         while(out(i) <= -pi)
%             out(i) = out(i) + 2*pi;
%         end
%     end
end

