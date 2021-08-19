function [] = makeLaserSaturationPlot()
%% Comments
% This script will plot the FWHM of the measure I II lineshape vs. Laser
% intensity for the data in Matt's capstone.
%
% T.E. Steinebrger
% August 2021

%% Script

%----% USER EDIT %----%
FILEDIRECTORY                                               = 'Results\';
SAVEDIRECTORY                                               = 'Results\';
POWERS                                                      = [0.0022 0.0022 0.0017 0.001 0.0007 0.0005]; % Watts
FILENUMBERS                                                 = [112 113 124 119 120 122]; 
%----% END USER EDIT %----%

% check for save directory
if ~exist(SAVEDIRECTORY, 'dir')
    mkdir(SAVEDIRECTORY);
end

% make laser intensity array
INTENSITIES                                                 = (2.*POWERS)/(pi*0.00025^2);
FWHM                                                        = zeros(1, length(INTENSITIES));
FILESTD                                                     = [0 0 0 0 0 0.1895]; % manually added by T.E.S.
COUNTER                                                     = 1;

for i = FILENUMBERS
    FWHM(COUNTER)                       = findFWHM(i, FILEDIRECTORY);
    COUNTER                                                 = COUNTER + 1;
end
ERROR                                                       = sqrt(((0.120*ones(1, length(FWHM))./FWHM).^2) + ((FILESTD./FWHM).^2));

figure('Position', [50 50 800 600], 'Color', 'white')
box on
semilogx(INTENSITIES, FWHM, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1.5)
hold on
errorbar(INTENSITIES, FWHM, ERROR.*FWHM, ERROR.*FWHM, 'LineStyle', 'none', 'Color', 'black', 'LineWidth', 1.5)
ylim([0 4.5])
xlim([3e3 3e4])
xlabel('I_{Las.} [W/m^2]', 'FontSize', 18)
ylabel('Spectral Width [GHz]', 'FontSize', 18)
set(gca, 'FontSize', 18, 'LineWidth', 2, 'ytick', [0 1 2 3 4])

% save directory
print([SAVEDIRECTORY 'LaserSaturation'], '-dpng')

% close
close all
end

%% FWHM function
function [FWHM] = findFWHM( varargin )
%% Comments
% this script will take the data from Matt's capstone and determine the 
% 'FWHM'.
% 
% T.E. Steinberger
% August 2021

%% Script
% check varargin
if isempty( varargin )
    FILENUMBER                                              = input('Enter File Number:\n');
    FILEDIRECTORY                                           = rawinput('EnterFileDirectory:\n');
else
    FILENUMBER                                              = varargin{1};
    FILEDIRECTORY                                           = varargin{2};
end

FILE                                                        = dir([FILEDIRECTORY '*' num2str(FILENUMBER) '\*' num2str(FILENUMBER) '.mat']);
load([FILE.folder '\' FILE.name], 'NORMALIZEDSIGNAL', 'RELATIVEFREQUENCY')

% find the indecies for FWHM
[~, MAXINDEX]                                               = max(NORMALIZEDSIGNAL);
LEFT                                                        = MAXINDEX - 550;
RIGHT                                                       = MAXINDEX + 200;
if LEFT < 0
    LEFT                                                    = 1;
elseif RIGHT > length(RELATIVEFREQUENCY)
    RIGHT                                                   = length(RELATIVEFREQUENCY);
end

[~, INDEX1]                                                 = min(abs(0.5 - NORMALIZEDSIGNAL(LEFT:(LEFT+200))));
[~, INDEX2]                                                 = min(abs(0.5 - NORMALIZEDSIGNAL(MAXINDEX:RIGHT)));
if FILENUMBER == 122
    FWHM                                                    = abs(RELATIVEFREQUENCY(INDEX1 + LEFT) - 0.4012); % hard code by T.E.S.
else
    FWHM                                                    = abs(RELATIVEFREQUENCY(INDEX1 + LEFT) - RELATIVEFREQUENCY(INDEX2 + MAXINDEX));
end

figure('Position', [50 50 800 600], 'Color', 'white')
box on
hold on
plot(RELATIVEFREQUENCY, NORMALIZEDSIGNAL, 'LineWidth', 2, 'Color', 'black')
plot([RELATIVEFREQUENCY(INDEX1 + LEFT); RELATIVEFREQUENCY(INDEX2 + MAXINDEX)], [NORMALIZEDSIGNAL(INDEX1 + LEFT); NORMALIZEDSIGNAL(INDEX2 + MAXINDEX)], 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red')
xlabel('\Deltaf [GHz]', 'FontSize', 18)
ylabel('Signal [arb.]', 'FontSize', 18)
set(gca, 'FontSize', 18, 'LineWidth', 1.5)

% save 
print(['Results\FWHM-' num2str(FILENUMBER)], '-dpng')

% close
close all
end