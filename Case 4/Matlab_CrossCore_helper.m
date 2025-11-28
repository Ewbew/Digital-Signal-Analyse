
% Eksempel på skrivning af array til fil, som kan læses ind i CrossCore
signal= [0.13 0.25 0.37 0.48 0.50 0.00 1.00 -1.00]; % værdier skal ligge mellem -1.00 og 1.00
powShort = 2^15; % Samme som shift med 15 bits
fid=fopen('minFil.dat', 'w');
for i=1:length(signal)-1
    fprintf(fid, '%d,\n', toRealnumber(signal(i), powShort)); 
end
fprintf(fid, '%d\n', toRealnumber(signal(i+1), powShort));   % Uden ","-tegn
fclose(fid);

% Funktion for convertering fra double/float til heltal
function y = toRealnumber(x, multFactor)
    y = round(x*multFactor);
    if (y >= multFactor) 
        y = multFactor-1; % Limitering af max værdi
    end
end




