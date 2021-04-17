% function [planes, offset] = my_sino_config(nring)
function [planes] = my_ring_config_uih(nring)
% my sinogram ordering 

% projection size per sinogram
% offset = zeros(nring, nring);
planes = zeros(nring*nring, 2);

c = 1;
for n = 1 : nring
	
	if n == 1

		for shift = 0 : (nring-n)
			% shift image
			planes(c,:) = [1+shift, n+shift];
			% offset(n+shift, 1+shift) = c;
			%%planes(c,:) = [1+shift, n+shift];
			%%offset(1+shift, n+shift) = c;
			c = c + 1;
			
		end

	else

		for shift = 0 : (nring-n)
			% shift image
			planes(c,:) = [1+shift, n+shift];
			% offset(n+shift, 1+shift) = c;
			%%planes(c,:) = [1+shift, n+shift];
			%%offset(1+shift, n+shift) = c;
			c = c + 1;
			
		end

		for shift = 0 : (nring-n)
			% shift image
			planes(c,:) = [n+shift, 1+shift];
			% offset(n+shift, 1+shift) = c;
			%%planes(c,:) = [1+shift, n+shift];
			%%offset(1+shift, n+shift) = c;
			c = c + 1;
			
			% if n > 1

			% 	planes(c,:) = [nring-n-shift+1, nring-shift];
			% 	% offset(nring-n-shift+1, nring-shift) = c;
			% 	%%planes(c,:) = [nring-shift, nring-n-shift+1];
			% 	%%offset(nring-shift, nring-n-shift+1) = c;
			% 	c = c + 1;
		
			% end
		end

	end



end

