

Sequence_length = 30;

% bit = randi([0 1], 1, Sequence_length);
% bit = [1 1 0 1 1 0 0 0 1 0 0];
% bit = [1 0 1 0 0 1 1 1 0 0 1];
bit = [0 0 1 1 1 0 1 0 0 0 0 1 0 1 0 1 0 1 1 1 0 0 0 1 0 1 1 1 1 0];
len = length(bit);
t = 0 : len ;

% Polar Non-return to Zero
nrz = 2 * bit - 1; 
nrz = [nrz nrz(len)];

% Non-return to Zero Inverted
nrzi = zeros(1, length(t)); % 0 -> nochange, 1 -> toggle
nrzi_bit = 1;

% Polar Return to Zero
rz = 0;
t_rz = 0 : 0.5 : len;

% Alternative Mark inversion (AMI)
ami = zeros(1, length(t)); %  BiPlor Non-return to Zero
ami_bit = 1;

% Manchester
manchester = 0;
t_manchester = 0 : 0.5 : len;

% Multi-level Transmission 3
mlt3 = zeros(1, length(t));% 0 -> nochange, 1 -> 1 or 0 or -1 (respectively)
mlt3_bit = 1;

for i = 1 : len
    % Non-return to Zero Inverted
    if bit(i)
        nrzi_bit = -1 * nrzi_bit;
        nrzi (i) = nrzi_bit;
    else
        nrzi(i) = nrzi_bit;
    end
    
    % Polar Return to Zero
    if bit(i)
        rz(i*2-1:i*2) = [1, 0];
    else
        rz(i*2-1:i*2) = [-1, 0];
    end
    
    % Alternative Mark inversion (AMI)
    if bit(i)
        ami(i) = ami_bit;
        ami_bit = -1 * ami_bit;
    end
    
    % Manchester
    if bit(i)
        manchester(i*2-1:i*2) = [1, -1];
    else
        manchester(i*2-1:i*2) = [-1, 1];
    end
    
    % Multi-level Transmission 3
    if bit(i)
        if mlt3_bit == 1
            mlt3(i) = 1;
            mlt3_bit = 2;
        elseif mlt3_bit == -1
            mlt3(i) = -1;
            mlt3_bit = 3;
        else
            mlt3(i) = 0;
            if mlt3_bit == 2
                mlt3_bit = -1;
            else
                mlt3_bit = 1;
            end
        end
    else
        if i ~= 1
            mlt3(i) = mlt3(i - 1);
        end
    end
end
rz = [rz 0];
manchester = [manchester 0];
        
% Power Spectral Density
fs = 100; % Sampling Frequency

[S_nrz, f_nrz] = pwelch(nrz, [], [], [], fs);
[S_nrzi, f_nrzi] = pwelch(nrzi, [], [], [], fs);
[S_rz, f_rz] = pwelch(rz, [], [], [], fs);
[S_ami, f_ami] = pwelch(ami, [], [], [], fs);
[S_manchester, f_manchester] = pwelch(manchester, [], [], [], fs);
[S_mlt3, f_mlt3] = pwelch(mlt3, [], [], [], fs);



% Plot

% Polar Non-return to Zero
subplot(6, 2, 1);
stairs(t, nrz, 'linewidth', 3);
title(num2str(bit));
set(gca, 'XTick', (0:length(t))); % adjust X-axis ticks
title('Polar Non-return to Zero');
grid on;
subplot(6, 2, 2);
plot(f_nrz, 10*log10(S_nrz), 'linewidth', 2);
title('PSD of Polar Non-return to Zero');
% ylim([-inf inf]);

% Non-return to Zero Inverted
subplot(6, 2, 3);
stairs(t, nrzi, 'linewidth', 3);
set(gca, 'XTick', (0:length(t)));
title('Non-return to Zero Inverted');
grid on;
subplot(6, 2, 4);
plot(f_nrzi, 10*log10(S_nrzi), 'linewidth', 2);
title('PSD of Non-return to Zero Inverted');
% ylim([-inf inf]);

% Polar Return to Zero
subplot(6, 2, 5);
stairs(t_rz, rz, 'linewidth', 3);
set(gca, 'XTick', (0:length(t)));
title('Polar Return to Zero');
grid on;
subplot(6, 2, 6);
plot(f_rz, 10*log10(S_rz), 'linewidth', 2);
title('PSD of Polar Return to Zero');
ylabel('Power (dB)');
% ylim([-inf inf]);

% Alternative Mark inversion (AMI)
subplot(6, 2, 7);
stairs(t, ami, 'linewidth', 3);
set(gca, 'XTick', (0:length(t)));
title('Alternative Mark inversion (AMI)');
grid on;
subplot(6, 2, 8);
plot(f_ami, 10*log10(S_ami), 'linewidth', 2);
title('PSD of Alternative Mark inversion (AMI)');
% ylim([-inf inf]);

% Manchester
subplot(6, 2, 9);
stairs(t_manchester, manchester, 'linewidth', 3);
set(gca, 'XTick', (0:length(t)));
title('Manchester');
grid on;
subplot(6, 2, 10);
plot(f_manchester, 10*log10(S_manchester), 'linewidth', 2);
title('PSD of Manchester');
% ylim([-inf inf]);

% Multi-level Transmission 3
subplot(6, 2, 11);
stairs(t, mlt3, 'linewidth', 3);
set(gca, 'XTick', (0:length(t)));
title('Multi-level Transmission 3');
grid on;
subplot(6, 2, 12);
plot(f_mlt3, 10*log10(S_mlt3), 'linewidth', 2);
title('PSD of Multi-level Transmission 3');
xlabel('Frequency (Hz)');
% ylim([-inf inf]);

set(gcf, 'Position', get(0, 'Screensize'));
suptitle(['Sequence: [', num2str(bit), ']']);