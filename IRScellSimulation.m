% This sets up and performs a simulation of a specific type of retinal
% ganglion cell (observed in mouse retina) under saccade-like image transitions.
% The simulation accompanies the manuscript "Sensitivity to image
% recurrence across eye-movement-like image transitions through local
% serial inhibition in the retina" by Krishnamoorthy, Weick, and Gollisch.
% Code written by T. Gollisch. This version from October 2016.

clear;
close all;

% some flags that can be used to run the simulation in different modes
saccadeMode = true;  % actual saccade transitions (true) or gray-screen transitions (false)
noBlack = false;     % only white bars are used when set to true
noWhite = false;     % only black bars are used when set to true
justFlash = false;   % stimulus during first fixation and during transition is set to gray when true
increasedTonicInhibition = false;  % increases the tonic inhibition of the On-type amacrine cell by a factor of 4 when set to true

% define time step and time duration (all in ms) for simulation
% 0-800 ms is first fixation (long enough so that transients from starting
% the simulation decay to zero before transition)
% 800-900 ms is the transition
% 900-1700 ms is the second fixation
dt = 0.1;
t = 0:dt:1700;
Nt = length(t);

% can be used to suppress output figures
figFilter = figure('Visible','on');
figStim = figure('Visible','on');
figOffBip = figure('Visible','on');
figOnBip = figure('Visible','on');
figOffAm = figure('Visible','on');
figOnAm = figure('Visible','on');
figInputs = figure('Visible','on');
figOffganglion = figure('Visible','on');

colors = ['c', 'b', 'g', 'r'];
set(0,'DefaultAxesTickDir', 'out');

% set filter for bipolar cells
% here for Off-type cell; On-type cell is the same with opposite sign
% (implemented during actual simulation)
Tfilterduration = 500;
Tfilter = 400;
Tdelay = 0;
tfilter = 0:dt:Tfilterduration;
filter = zeros(size(tfilter));
filter_part1 = zeros(size(tfilter));
filter_part2 = zeros(size(tfilter));
filter_part1(tfilter>Tdelay) = exp(-(tfilter(tfilter>Tdelay)-Tdelay).^2/(0.1*(Tfilter-Tdelay))^2)-exp(-(tfilter(tfilter>Tdelay)-Tdelay).^2/(0.2*(Tfilter-Tdelay))^2);
filter_part2(tfilter>Tdelay) = 0.5*( exp(-(tfilter(tfilter>Tdelay)-Tdelay).^2/(0.2*(Tfilter-Tdelay))^2)-exp(-(tfilter(tfilter>Tdelay)-Tdelay).^2/(0.4*(Tfilter-Tdelay))^2) );
filter_part2 = filter_part2*(sum(filter_part1)/sum(filter_part2)); % this makes sure that the positive and negative lobes of the filter cancel out so that the activation is transient
filter = filter_part1-filter_part2;
normfactor = sqrt(sum(filter.*filter));
filter = filter/normfactor;

% set low-pass filter for On-type amacrine cell
Tlowpass = 600;
Tdelaylowpass = 0;
tlowpass = 0:dt:Tlowpass;
lowpass = zeros(size(tlowpass));
lowpass(tlowpass>Tdelaylowpass) = -exp(-(tlowpass(tlowpass>Tdelaylowpass)-Tdelaylowpass).^2/(0.05*Tlowpass)^2)+exp(-(tlowpass(tlowpass>Tdelaylowpass)-Tdelaylowpass).^2/(0.15*Tlowpass)^2);
normfactor = sqrt(sum(lowpass.*lowpass));
lowpass = lowpass/normfactor;

% plot the filters
set(0,'CurrentFigure',figFilter);
plot( tfilter, filter, 'b' );
hold on;
plot( tlowpass, lowpass, 'r' );
hold off;

% define the stimulus; bright = 1; dark = -1
% Xstim and Ystim are grating contrasts in the four subunits;
% Xstim for first fixation; Ystim for second fixation
Xstim = [1; 1; -1; -1];
Ystim = [1; 1; -1; -1];
Xstim = circshift( Xstim, -1 ); % this is just to undo the first circular permution at the beginning of the for-loops below
Ystim = circshift( Ystim, -1 );
if(noBlack)
    Xstim(find(Xstim<0)) = 0;
    Ystim(find(Ystim<0)) = 0;
end;
if(noWhite)
    Xstim(find(Xstim>0)) = 0;
    Ystim(find(Ystim>0)) = 0;
end;

for x=1:4
    Xstim = circshift( Xstim, 1 ); % do circular permutation of grating contrasts, corresponding to a quarter-period shift
    for y=1:4
        Ystim = circshift( Ystim, 1 );
        
        % the variables that hold the stimulus and neuronal activations
        % over the stimulus duration
        stimulus = zeros( 4, Nt );
        OFFbipolar = zeros( 4, Nt );
        ONbipolar = zeros( 4, Nt );
        OFFamacrine = zeros( 4, Nt );
        ONamacrine = zeros( 4, Nt );
        OFFganglionInput = zeros( 4, Nt );
        OFFganglion = zeros( 1, Nt );
        
        % set stimulus
        for i=1:4  % looping through the four subunits
            stimulus( i, t<=800 ) = Xstim(i);
            stimulus( i, t>900 ) = Ystim(i);
            if( saccadeMode ) % in saccade mode, the stimulus for each subunit during the transition oscillates in a triangular fashion between 1 and -1
                Nsubdivisions = 8 + (y-x);
                if( Nsubdivisions < 6 )
                    Nsubdivisions = Nsubdivisions + 4;
                end;
                if( Nsubdivisions > 10 )
                    Nsubdivisions = Nsubdivisions - 4;
                end;
                currentIndex = i;
                newIndex = i-1;
                if( newIndex < 1 )
                    newIndex = newIndex + 4;
                end;
                for j=1:Nsubdivisions
                    steps = 1:length( stimulus( i, (t>800+(j-1)*100.0/Nsubdivisions) & (t<=800+j*100.0/Nsubdivisions) ) );
                    steps = (steps-1.0)/length(steps);
                    stimulus( i, (t>800+(j-1)*100.0/Nsubdivisions) & (t<=800+j*100.0/Nsubdivisions) ) = Xstim(currentIndex)*steps(end:-1:1)+Xstim(newIndex)*steps;
                    currentIndex = newIndex;
                    newIndex = newIndex - 1;
                    if( newIndex<1 )
                        newIndex = newIndex + 4;
                    end;
                end;
            else % for a masked saccade, the stimulus is set to zero during the transition
                stimulus(i, (t>800) & (t<=900) ) = 0;
            end;
            
            if( justFlash ) % if used, this sets all stimulus values to zero (gray) before second fixation
                stimulus(i, t<=900 ) = 0;
            end;
        end;
        
        for i=1:4 % loop through all four subunits
            signal = conv( filter, stimulus(i,:) );
            OFFbipolar(i,:) = signal(1:Nt);
            ONbipolar(i,:) = -1*OFFbipolar(i,:);
            %ONbipolar(i,:) = 0;  % simulating L-AP4 when uncommented
            OFFamacrine(i,:) = OFFbipolar(i,:);
            OFFamacrine(i,find(OFFamacrine(i,:)<0)) = 0; % rectificatio of input into the Off-type amacrine cell
            %OFFamacrine(i,:) = 0; % simulating gabazine when uncommented
            ONamacrine(i,:) = 0.025*ONbipolar(i,:);
            ONamacrine(i,find(ONamacrine(i,:)<0)) = 0; % rectification of input into the On-type amacrine cell
            ONamacrine(i,:) = ONamacrine(i,:)-0.1*OFFamacrine(i,:); % inhibition from the Off-type to the On-type amacrine cell; note that the Off-type amacrine cell never has negative activation, so no rectification necessary in implementation
            %low pass filter of ONamacrine
            signal = conv( lowpass, ONamacrine(i,:) );
            ONamacrine(i,:) = signal(1:Nt); % cuts the activation vector back to the original length
            tonicInhibition = 2;
            if( increasedTonicInhibition )
                tonicInhibition = tonicInhibition * 4;
            end;
            ONamacrine(i,:) = ONamacrine(i,:)+tonicInhibition;  % tonic inhibition implemented by sustained activity of On-type amacrine cell
            %ONamacrine(i,:) = 0; % simulating strychnine when uncommented
            ONamacrine(i,find(ONamacrine(i,:)<0)) = 0; % rectification of the signal transmitted by the On-type amacrine cell
        end;
        for i=1:4 % loop through all four subunits to sum the local inputs to the ganglion cell
            OFFganglionInput(i,:) = OFFbipolar(i,:);
            for j=1:4 % presynaptic inhibition from the On-type amacrine cells onto the Off-bipolar cell terminal
                OFFganglionInput(i,:) = OFFganglionInput(i,:) - 0.5*ONamacrine(j,:);
            end;
            OFFganglionInput(i,find(OFFganglionInput(i,:)<0)) = 0;  % rectification of the excitatory input
            OFFganglionInput(i,:) = OFFganglionInput(i,:) - 1*ONamacrine(i,:); % postsynaptic inhibition that acts as an independent input onto the ganglion cell
            OFFganglion = OFFganglion + OFFganglionInput(i,:); % summation by the ganglion cell of its inputs
        end;
        OFFganglion(find(OFFganglion<0))=0; % rectification of ganglion cell activity
        
        % a few plots for monitoring the stimulus and the activity of
        % different neurons in the four subunits for the 16 transitions:
        
        % stimulus
        set(0,'CurrentFigure',figStim);
        subplot(4,4,(4-x)*4+y);
        for i=1:4
            hold on;
            plot( t, stimulus(i,:)+(i-1)*3, colors(i) );
            hold off;
        end;

        % Off-type bipolar cells
        set(0,'CurrentFigure',figOffBip);
        subplot(4,4,(4-x)*4+y);
        for i=1:4
            hold on;
            plot( t, max(OFFbipolar(i,:),0)+(i-1)*30, colors(i) );
            hold off;
        end;

        % On-type bipolar cells
        set(0,'CurrentFigure',figOnBip);
        subplot(4,4,(4-x)*4+y);
        for i=1:4
            hold on;
            plot( t, ONbipolar(i,:)+(i-1)*30, colors(i) );
            hold off;
        end;

        % Off-type amacrine cells
        set(0,'CurrentFigure',figOffAm);
        subplot(4,4,(4-x)*4+y);
        for i=1:4
            hold on;
            plot( t, OFFamacrine(i,:)+(i-1)*30, colors(i) );
            hold off;
        end;

        % On-type amacrine cells
        set(0,'CurrentFigure',figOnAm);
        subplot(4,4,(4-x)*4+y);
        for i=1:4
            hold on;
            plot( t, ONamacrine(i,:)+(i-1)*30, colors(i) );
            hold off;
        end;

        % summed local inputs to the ganglion cell
        set(0,'CurrentFigure',figInputs);
        subplot(4,4,(4-x)*4+y);
        for i=1:4
            hold on;
            plot( t, OFFganglionInput(i,:)+(i-1)*30, colors(i) );
            hold off;
        end;

        % reduce sampling of curves for easier plotting of the actual model
        % output
        OFFganglion_plotting = decimate(OFFganglion, 20);
        t_plotting = decimate(t, 20);
        
        % firing rate of ganglion cell
        set(0,'CurrentFigure',figOffganglion);
        subplot(4,4,(4-x)*4+y);
        plot( t_plotting, OFFganglion_plotting(:), 'k', 'LineWidth', 2 );
        hold on
        plot([800 800],[-2 37.5],'k','LineWidth',1);
        plot([900 900],[-2 37.5],'k','LineWidth',1);
        axis([800 1200 -2 37.5]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        axis off;

    end;
end;
