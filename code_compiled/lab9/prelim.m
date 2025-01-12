%% Prelim Q2
% Determine all possible combinations of a, p, q to get desired a'

% test function
combinations(1);

%% create image
im1 = fourier(128,64,32,32,50,250,200);
im2 = fourier(50,50,50,50,50,50,50);
im3 = fourier(50,50,50,50,20,20,20);
im4 = fourier(0,255,0,0,30,1,1);
figure
subplot(3,4,1)
imagesc(im1);
subplot(3,4,2)
imagesc(im2);
subplot(3,4,3)
imagesc(im3);
subplot(3,4,4)
imagesc(im4);
set(gcf,'color','w');


% fourier transform of the image
IfftTemp = fft2(im1); 
it2 = fft(im2);
it3 = fft(im3);
it4 = fft(im4);

% Take abs val and undo the shift induced by fft
Ifft = fftshift(abs(IfftTemp));
ift2 = fftshift(abs(it2));
ift3 = fftshift(abs(it3));
ift4 = fftshift(abs(it4));

subplot(3,4,5)
imagesc(Ifft, [1e-13 1e-9]);
subplot(3,4,6)
imagesc(ift2, [1e-13 1e-9]);
subplot(3,4,7)
imagesc(ift3, [1e-13 1e-9]);
subplot(3,4,8)
imagesc(ift4, [1e-13 1e-9]);

% experimental ones


%% FUNCTION: takes in the desired value for a', 
% evaluates all possible combinations for p, q, and a and returns 
% an array with rows a, p , q: in centimeters
function choose = combinations(aprime)
    % wavelength in centimeters
    lambda = 6.25e-5;
    % import files with a and p
    a = xlsread('Possible_a_Values.xlsx');
    a = a./1e4; % convert from microm. to cm
    a = a.';
    q = xlsread('Possible_q_Values.xlsx');
    q = q.';
    q = q./10; % convert from mm to cm
    
    % take all possible combinations of the arrays
    combs = combvec(a,q);
    % store combinations that work
    choose = zeros(3,1);
    
    % take product of all vectors
    colcount = 0;
    for i=1:length(combs)
        a0 = combs(1,i);
        q0 = combs(2,i);
        
        % p will be in cm
        p0 = ((aprime*sqrt(lambda)/a0)-(1/q0))^-1;
        
        % if p is within range
        if (p0 > 1.5 && p0 < 150)
            colcount = colcount + 1;
            choose(:,colcount) = [a0;p0;q0];
        end
    end
end

%% FUNCTION: fourier transform

function im = fourier(a0, a1, a2, a3, l1, l2 , l3)
    % array to store
    im = zeros(600,800);
    for x = 1:600
        for y = 1:800
            im(x,y) = a0 + a1*sin((2*pi*y)/l1) + a2*sin((2*pi*y)/l2) + ...
                a3*sin((2*pi*x)/l3);
        end
    end
end