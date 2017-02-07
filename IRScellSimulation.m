% This sets up and performs a simulation of a specific type of retinal
% ganglion cell (observed in mouse retina) under saccade-like image transitions.
% The simulation accompanies the manuscript "Sensitivity to image
% recurrence across eye-movement-like image transitions through local
% serial inhibition in the retina" by Krishnamoorthy, Weick, and Gollisch.
% Code written by T. Gollisch. This version from February 2017.

clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% setting flags to control simulation mode %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This can be used to run the simulation in different modes as used in the
% manuscript; default is "false" for all flags
maskedSaccade = false;  % actual saccade transitions (false) or gray-screen transitions (true)
noBlack = false;     % only white bars are used when set to true
noWhite = false;     % only black bars are used when set to true
justFlash = false;   % stimulus during first fixation and during transition is set to gray when true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% setting simulation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tTotal = 1700.0;    % default: 1700; total simulated time (in ms)
tDuration = 100.0;  % default: 100; duration of the transition
tOnset = 800.0;     % default: 800; time when the transition starts

dt = 0.1;         % default: 0.1; time step size used in simulation
tPlotStart = 800; % default: 800; setting the start of the plot windows; for viewing traces from interneurons, default: 700
tPlotEnd = 1200;  % default: 1200; setting the end of the plot windows; for viewing traces from interneurons, default: 1500

Tfilterduration = 600; % default: 600; duration over which filters are considered; should be long enough to capture all relevant parts of the filters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% free parameters of the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T1 = 40;   % default: 40; time constant of excitatory lobe of bipolar cell filter
T2 = 60;   % default: 60; time constant of suppressive lobe of bipolar cell filter

Tlowpass_rise = 30;   % default: 30; first time constant of low-pass filter of On-type amacrine cell
Tlowpass_decay = 200; % default: 200; second time constant of low-pass filter

weight_ONBC_ONAC = 0.05;  % default: 0.05; synaptic weight of the connection from On-type BC to On-type AC
weight_OFFAC_ONAC = 0.2;  % default: 0.2; synaptic weight of the connection from Off-type AC to On-type AC

tonicInhibition = 2.0;   % default 2.0; level of tonic activation of ON-type AC; default for simulation with increased tonic inhibition: 16.0

relativeStrengthOfPresynapticInhibition = 0;
% the On-type inhibition onto the GC can be postsynaptic or presynaptic,
% the relative presynaptic strength can be set here
% should be a number between 0 (no presynaptic inh.) and 1 (only presyn.)
% GC firing rates should not be affected, but plots of local inputs will be


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% defining the plot windows %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figFilter = figure('Visible','on');
figStim = figure('Visible','on');
figOffBip = figure('Visible','on');
figOnBip = figure('Visible','on');
figOffAm = figure('Visible','on');
figOnAm = figure('Visible','on');
figInputs = figure('Visible','on');
figSampleTraces = figure('Visible','on');
figOffganglion = figure('Visible','on');
% parameters can be used to suppress output figures

colors = ['c', 'b', 'g', 'r']; % setting the order of used colors for the traces in the four subunits
set(0,'DefaultAxesTickDir', 'out');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% defining the filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set filter for bipolar cells
% here for Off-type cell; On-type cell is the same with opposite sign
% (implemented during actual simulation)
tfilter = 0:dt:Tfilterduration;
filter_part1 = exp(-tfilter.^2/T1^2)-exp(-tfilter.^2/(2*T1)^2);
filter_part2 = exp(-tfilter.^2/T2^2)-exp(-tfilter.^2/(2*T2)^2);
filter_part2 = filter_part2*(sum(filter_part1)/sum(filter_part2)); % this makes sure that the positive and negative lobes of the filter cancel out so that the activation is transient
filter = filter_part1-filter_part2;
normfactor = sqrt(sum(filter.*filter));
filter = filter/normfactor;

% set low-pass filter for On-type amacrine cell
tlowpass = 0:dt:Tfilterduration;
lowpass = exp(-tlowpass.^2/Tlowpass_decay^2)-exp(-tlowpass.^2/Tlowpass_rise^2);
normfactor = sqrt(sum(lowpass.*lowpass));
lowpass = lowpass/normfactor;

% plot the filters
set(0,'CurrentFigure',figFilter);
plot( tfilter, filter, 'b' );
hold on;
plot( tlowpass, lowpass, 'r' );
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% some preparations for the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the sampling time points for the simulation
t = 0:dt:tTotal;
Nt = length(t);

% define the stimulus; bright = 1; dark = -1
% Xstim and Ystim are grating contrasts in the four subunits;
% Xstim for starting image; Ystim for target image
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% loop through all 16 transitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for x=1:4   % 4 different starting positions
    Xstim = circshift( Xstim, 1 ); % do circular permutation of grating contrasts, corresponding to a quarter-period shift
    for y=1:4   % 4 different target positions
        Ystim = circshift( Ystim, 1 );
        
        % the variables that hold the stimulus and neuronal activations
        % over the stimulus duration
        stimulus = zeros( 4, Nt );
        OFFbipolar = zeros( 4, Nt );
        ONbipolar = zeros( 4, Nt );
        OFFamacrine = zeros( 4, Nt );
        ONamacrine = zeros( 4, Nt );
        ONamacrineSynapse = zeros( 4, Nt );
        OFFganglionInput = zeros( 4, Nt );
        OFFganglion = zeros( 1, Nt );
        
        % set stimulus
        for i=1:4  % looping through the four subunits
            stimulus( i, t <= tOnset ) = Xstim(i);
            stimulus( i, t > tOnset+tDuration ) = Ystim(i);
            if( maskedSaccade )
                % for a masked saccade, the stimulus is set to zero during the transition
                stimulus(i, (t>tOnset) & (t<=tOnset+tDuration) ) = 0;
            else
                % in saccade mode, the transition is partitioned into subdivisions
                % each full cycle corresponds to 4 subdivisions
                % in each subdivision, one stimulus subunit goes linearly
                % from bright to dark, one goes linearly from dark to
                % bright, and two stay constant at dark and bright,
                % respectively
                Nsubdivisions = 8 + (y-x);
                % this corresponds to two full cylces (8 subdivisions), plus corrections to account for the difference between starting and target position
                if( Nsubdivisions < 6 )
                    Nsubdivisions = Nsubdivisions + 4; % this ensures that the transition goes through at least 1.5 cycles (6 subdivisions)
                end;
                if( Nsubdivisions > 10 )
                    Nsubdivisions = Nsubdivisions - 4; % this ensures that the transition goes through at most 2.5 cycles (10 subdivisions)
                end;
                % typically, we use around 2 full cycles of the grating
                % during the transition, but for short transitions (here
                % less than 50 ms), we use only around 1 cycle, so the
                % following adjustment takes place
                if( tDuration<50 )
                    Nsubdivisions = 4 + (y-x);
                    if( Nsubdivisions < 2 )
                        Nsubdivisions = Nsubdivisions + 4;
                    end;
                    if( Nsubdivisions > 6 )
                        Nsubdivisions = Nsubdivisions - 4;
                    end;
                end;
                currentIndex = i; % index for the stimulus value at the beginning of the subdivision
                newIndex = i-1;   % index for the stimulus value at the end of the subdivision
                if( newIndex < 1 )
                    newIndex = newIndex + 4;
                end;
                for j=1:Nsubdivisions
                    steps = 1:length( stimulus( i, (t>tOnset+(j-1)*tDuration/Nsubdivisions) & (t<=tOnset+j*tDuration/Nsubdivisions) ) );
                    steps = (steps-1.0)/length(steps);
                    stimulus( i, (t>tOnset+(j-1)*tDuration/Nsubdivisions) & (t<=tOnset+j*tDuration/Nsubdivisions) ) = Xstim(currentIndex)*steps(end:-1:1)+Xstim(newIndex)*steps;
                    currentIndex = newIndex;
                    newIndex = newIndex - 1;
                    if( newIndex<1 )
                        newIndex = newIndex + 4;
                    end;
                end;
            end;
            
            if( justFlash ) % if used, this sets all stimulus values to zero (gray) before second fixation
                stimulus(i, t<=tOnset+tDuration ) = 0;
            end;
        end;
        
        % do the actual simulation for all cells in each subunit
        for i=1:4 % loop through all four subunits
            signal = conv( filter, stimulus(i,:) );
            
            OFFbipolar(i,:) = signal(1:Nt);
            ONbipolar(i,:) = -1*OFFbipolar(i,:);
            
            OFFamacrine(i,:) = OFFbipolar(i,:);
            OFFamacrine(i,find(OFFamacrine(i,:)<0)) = 0; % rectification of input into the Off-type amacrine cell
            %OFFamacrine(i,:) = 0; % simulating gabazine when uncommented
            
            ONamacrine(i,:) = weight_ONBC_ONAC * ONbipolar(i,:);
            ONamacrine(i,find(ONamacrine(i,:)<0)) = 0; % rectification of input into the On-type amacrine cell

            ONamacrine(i,:) = ONamacrine(i,:) - weight_OFFAC_ONAC * OFFamacrine(i,:);
            % inhibition from the Off-type to the On-type amacrine cell; note that the Off-type amacrine cell never has negative activation, so no rectification necessary in implementation            
            
            %low pass filter of ONamacrine
            signal = conv( lowpass, ONamacrine(i,:) );
            ONamacrine(i,:) = signal(1:Nt); % cuts the activation vector back to the original length
            
            ONamacrine(i,:) = ONamacrine(i,:) + tonicInhibition;  % tonic inhibition implemented by sustained activity of On-type amacrine cell
            
            ONamacrineSynapse(i,:) = ONamacrine(i,:);  % takes the output signal of the On-type AC to implement rectification
            ONamacrineSynapse(i,find(ONamacrineSynapse(i,:)<0)) = 0; % rectification of the signal transmitted by the On-type amacrine cell            
            %ONamacrineSynapse(i,:) = 0; % simulating strychnine when uncommented
        end;

        for i=1:4 % loop through all four subunits to sum the local inputs to the ganglion cell
            OFFganglionInput(i,:) = OFFbipolar(i,:);
            for j=1:4 % presynaptic inhibition from the On-type amacrine cells onto the Off-bipolar cell terminal
                OFFganglionInput(i,:) = OFFganglionInput(i,:) - relativeStrengthOfPresynapticInhibition*0.5*ONamacrineSynapse(j,:);
            end;
            OFFganglionInput(i,find(OFFganglionInput(i,:)<0)) = 0;  % rectification of the excitatory input
            OFFganglionInput(i,:) = OFFganglionInput(i,:) - (1-relativeStrengthOfPresynapticInhibition)*ONamacrineSynapse(i,:); % postsynaptic inhibition that acts as an independent input onto the ganglion cell
            OFFganglion = OFFganglion + OFFganglionInput(i,:); % summation by the ganglion cell of its inputs
        end;
        OFFganglion(find(OFFganglion<0)) = 0; % rectification of ganglion cell activity
        
        % a few plots for monitoring the stimulus and the activity of
        % different neurons in the four subunits for the 16 transitions:
        
        % stimulus
        set(0,'CurrentFigure',figStim);
        subplot(4,4,(4-x)*4+y);
        for i=1:4
            hold on;
            plot( t, stimulus(i,:)+(i-1)*3, colors(i) );
            hold off;
            xlim([tPlotStart tPlotEnd]);
        end;
        if(x==4 & y==1) title('Stimulus'); else title(''); end;

        % Off-type bipolar cells
        set(0,'CurrentFigure',figOffBip);
        subplot(4,4,(4-x)*4+y);
        for i=1:4
            hold on;
            plot( t, OFFbipolar(i,:)+(i-1)*30, colors(i) );
            hold off;
            xlim([tPlotStart tPlotEnd]);
        end;
        if(x==4 & y==1) title('Off BC'); else title(''); end;

        % On-type bipolar cells
        set(0,'CurrentFigure',figOnBip);
        subplot(4,4,(4-x)*4+y);
        for i=1:4
            hold on;
            plot( t, ONbipolar(i,:)+(i-1)*30, colors(i) );
            hold off;
            xlim([tPlotStart tPlotEnd]);
        end;
        if(x==4 & y==1) title('On BC'); else title(''); end;

        % Off-type amacrine cells
        set(0,'CurrentFigure',figOffAm);
        subplot(4,4,(4-x)*4+y);
        for i=1:4
            hold on;
            plot( t, OFFamacrine(i,:)+(i-1)*30, colors(i) );
            hold off;
            xlim([tPlotStart tPlotEnd]);
        end;
        if(x==4 & y==1) title('Off AC'); else title(''); end;

        % On-type amacrine cells
        set(0,'CurrentFigure',figOnAm);
        subplot(4,4,(4-x)*4+y);
        for i=1:4
            hold on;
            plot( t, ONamacrine(i,:)+(i-1)*30, colors(i) );
            hold off;
            xlim([tPlotStart tPlotEnd]);
        end;
        if(x==4 & y==1) title('On AC'); else title(''); end;

        % summed local inputs to the ganglion cell
        set(0,'CurrentFigure',figInputs);
        subplot(4,4,(4-x)*4+y);
        for i=1:4
            hold on;
            plot( t, OFFganglionInput(i,:)+(i-1)*30, colors(i) );
            hold off;
            xlim([tPlotStart tPlotEnd]);
        end;
        if(x==4 & y==1) title('Summed Local Inputs'); else title(''); end;
        
        % sample traces for four locations
        % with different dark/bright transitions
        if(x==1 & y==2)
            set(0,'CurrentFigure',figSampleTraces);
            for i=1:4
                subplot(2,2,i);
                hold on;
                % reduce sampling of curves for easier plotting
                OFFbipolar_plotting = decimate( OFFbipolar(i,:), 40 );
                ONamacrine_plotting = decimate( ONamacrine(i,:), 40 );
                OFFganglionInput_plotting = decimate( OFFganglionInput(i,:), 40 );
                t_plotting = decimate( t, 40 );
                plot( t_plotting, OFFbipolar_plotting, 'c', 'LineWidth', 3 );
                plot( t_plotting, ONamacrine_plotting, 'm', 'LineWidth', 3 );
                plot( t_plotting, OFFganglionInput_plotting, 'k', 'LineWidth', 2 );
                plot( [tPlotStart tPlotEnd], [0 0], 'k-' );
                plot([tOnset tOnset],[-50 50],'k-','LineWidth',1);
                plot([tOnset+tDuration tOnset+tDuration],[-50 50],'k-','LineWidth',1);
                hold off;
                xlim([tPlotStart tPlotEnd]);
                ylim([-50 50]);
                if(i==1) legend('OFF BC', 'ON AC', 'GC input'); end;
                switch i
                    case 1, title('White -- Transition -- Black');
                    case 2, title('White -- Transition -- White');
                    case 3, title('Black -- Transition -- White');
                    case 4, title('Black -- Transition -- Black');
                end;
                %set(gca,'xtick',[]);
                %set(gca,'ytick',[]);
                %axis off;
            end;
        end;
        
        % reduce sampling of curves for easier plotting
        OFFganglion_plotting = decimate(OFFganglion, 20);
        t_plotting = decimate(t, 20);
        
        % firing rate of ganglion cell
        set(0,'CurrentFigure',figOffganglion);
        subplot(4,4,(4-x)*4+y);
        plot( t_plotting, OFFganglion_plotting(:), 'k', 'LineWidth', 2 );
        hold on
        plot([tOnset tOnset],[-2.2 40],'k','LineWidth',1);
        plot([tOnset+tDuration tOnset+tDuration],[-2.2 40],'k','LineWidth',1);
        axis([tPlotStart tPlotEnd -2.2 40]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        axis off;

    end;
end;
