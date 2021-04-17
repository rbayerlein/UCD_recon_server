function make_plaeff

fname = '/mnt/data/rbayerlein/code/explorer-master/reconstruction/plaeff_451584.NC'; 

fid_p = fopen(fname, 'r'); 

plaeff = fread(fid_p, inf, 'float'); 


michel = zeros(672, 672); 

michel_lut = zeros(672, 672); 

m = 0;
x = 1; 
y = 1; 

counter = 1; 
for k = 1:672
	michel(k,k) = plaeff(counter);
	michel_lut(k,k) = counter-1; 
	counter = counter+1; 
end


mend = 671;  

x = 2; 
y = 1; 

for k = 2:(672)
	mend = 672 - k + 1; 
	x = k; 
	y = 1; 
	for kk = 1:mend

		michel(x,y) = plaeff(counter); 
		michel_lut(x,y) = counter-1;
		x = x + 1; 
		y = y + 1; 
		counter = counter + 1; 
	end

	y = k; 
	x = 1; 

	for kk = 1:mend
		michel(x,y) = plaeff(counter); 
		michel_lut(x,y) = counter-1;
		x = x + 1; 
		y = y + 1; 
		counter = counter + 1; 
	end


end



figure
imagesc(michel)


fname = 'michel_lut'; 
fid_o = fopen(fname, 'w'); 
fwrite(fid_o, michel_lut, 'int'); 

fclose(fid_o); 



