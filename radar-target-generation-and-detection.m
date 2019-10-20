clear all
close all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_range = 200;
d_res = 1;
max_vel = 100;
c = 3e8; %speed of light

%% User Defined Range and Velocity of target
% define the target's initial position and velocity. Note : Velocity
% remains contant
ini_pos = 100;
tvel = -40;


%% FMCW Waveform Generation

%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
% d_res = c / 2bw;
B = c / (2 * d_res);
Tchirp = (2 * max_range * 5.5) / c;
slope = B / Tchirp;

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq
 
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t) 
    %For each time stamp update the Range of the Target for constant velocity.
    tpos = ini_pos + (t(i) * tvel);
    delay = 2*tpos/c;

    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*((fc*t(i)) + (slope*t(i)^2)/2));
    Rx(i) = cos(2*pi*( fc*(t(i)-delay) + (slope*(t(i) - delay)^2)/2 ));

    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
end

%% RANGE MEASUREMENT

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix, [Nr, Nd]);
%run the FFT on the beat signal along the range bins dimension (Nr)
range_fft = fft(Mix,Nr);

% Take the absolute value of FFT output
range_fft = abs(range_fft);

% Normalize after taking the absolute value
range_fft = range_fft./max(range_fft);

% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
range_fft = range_fft(1:(Nr/2));

% plot FFT output 
figure ('Name','Range FFT')
plot(range_fft);
title("Range FFT");
axis ([0 200 0 1]);
ylabel('Amplitude (normalized)');
xlabel('Range (m)');


%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.
Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
%figure,surf(doppler_axis,range_axis,RDM);

figure ('Name','FFT2 Output - RDM 3D');
surf(doppler_axis,range_axis,RDM);
title('FFT2 Output - RDM 3D');
xlabel('Speed');
ylabel('Range');
zlabel('Amplitude');
colorbar;

% Additional plots
figure ('Name','Amplitude and Range');
surf(doppler_axis,range_axis,RDM);
title('FFT2 Range');
xlabel('Speed');
ylabel('Range');
zlabel('Amplitude');
view(90,0);

figure ('Name','Amplitude and Speed');
surf(doppler_axis,range_axis,RDM);
title('FFT2 Speed');
xlabel('Speed');
ylabel('Range');
zlabel('Amplitude');
view(0,0);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

%Select the number of Training Cells in both the dimensions.
Tr = 10;      %the number of Training Cells in both the dimensions.
Td = 8;

%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 4;       %the number of Guard Cells in both dimensions around the 
Gd = 4;       %Cell under test (CUT) for accurate estimation
 
% offset the threshold by SNR value in dB
offset = 10;

%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.

   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR
%RDM = RDM/max(max(RDM)); % Normalizing
 
training_size = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1) - (2*Gr+1)*(2*Gd+1);
CFAR2D = zeros(size(RDM));
threshM = zeros(size(RDM));
offsetM = zeros(size(RDM));
%disp(training_size);

% 2. Slide window across the signal length

for i = (1+Gr+Tr):((Nr/2)-(Gr+Tr))
    for j = (1+Gd+Td):((Nd)-(Gd+Td))
        %Create a vector to store noise_level for each iteration on training cells
        noise_level = zeros(1,1);
        % loop through the training cells
        for p = i-(Gr+Tr):i+(Gr+Tr)
            for q = j-(Gd+Td):j+(Gd+Td)
                %skip guard cells
                if(abs(i-p)>Gr || abs(j-q)>Gd)
                    noise_level = noise_level + db2pow(RDM(p,q));
                end
            end
        end
        % calculate average        
        threshold = pow2db(noise_level/training_size);
        threshM(i,j) = threshold;
        offsetM(i,j) = threshM(i,j) + offset;
        if (RDM(i,j) > threshold + offset)
            CFAR2D(i,j) = 1;
        else
            CFAR2D(i,j) = 0;
        end
    end
end

% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 

% not necessary, since CFAR2D matrix was initialized with zeroes

%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure ('Name','CFAR2D Output');
surf(doppler_axis,range_axis,CFAR2D);
title('CFAR Output');
xlabel('Speed');
ylabel('Range');
zlabel('Amplitude');
colorbar;

figure ('Name','CFAR Range');
surf(doppler_axis,range_axis,CFAR2D);
title('CFAR Range');
xlabel('Speed');
ylabel('Range');
zlabel('Amplitude');
view(90,0);

figure ('Name','CFAR Speed');
surf(doppler_axis,range_axis,CFAR2D);
title('CFAR Speed');
xlabel('Speed');
ylabel('Range');
zlabel('Amplitude');
view(0,0);


figure ('Name','Thresh');
surf(doppler_axis,range_axis,threshM);
title('Thresh');
xlabel('Speed');
ylabel('Range');
zlabel('Amplitude');
colorbar;

figure ('Name','Offset');
surf(doppler_axis,range_axis,offsetM);
title('Offset');
xlabel('Speed');
ylabel('Range');
zlabel('Amplitude');
colorbar;

disp((2*(Td+Gd+1)*2*(Tr+Gr+1)-(Gr*Gd)-1));
disp(training_size);
