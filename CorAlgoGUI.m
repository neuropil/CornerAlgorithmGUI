function varargout = CorAlgoGUI(varargin)
% CORALGOGUI MATLAB code for CorAlgoGUI.fig
%      CORALGOGUI, by itself, creates a new CORALGOGUI or raises the existing
%      singleton*.
%
%      H = CORALGOGUI returns the handle to a new CORALGOGUI or the handle to
%      the existing singleton*.
%
%      CORALGOGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORALGOGUI.M with the given input arguments.
%
%      CORALGOGUI('Property','Value',...) creates a new CORALGOGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CorAlgoGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CorAlgoGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CorAlgoGUI

% Last Modified by GUIDE v2.5 21-Mar-2014 16:03:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CorAlgoGUI_OpeningFcn, ...
    'gui_OutputFcn',  @CorAlgoGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before CorAlgoGUI is made visible.
function CorAlgoGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CorAlgoGUI (see VARARGIN)

% Choose default command line output for CorAlgoGUI
handles.output = hObject;

handles.image2show = imread('imageDraw.tif');
axes(handles.polyDraw)
imshow(handles.image2show);

handles.noSymbol = imread('NOsign.jpg');
handles.yesSymbol = imread('YESsign.png');

cla(handles.check_yTL)
cla(handles.check_xTL)
cla(handles.check_yBL)
cla(handles.check_xBL)
cla(handles.check_yBR)
cla(handles.check_xBR)
cla(handles.check_yTR)
cla(handles.check_xTR)

set(handles.check_yTL,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');
set(handles.check_xTL,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');
set(handles.check_yBL,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');
set(handles.check_xBL,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');
set(handles.check_yBR,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');
set(handles.check_xBR,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');
set(handles.check_yTR,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');
set(handles.check_xTR,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');



% TO check
% axes(handles.check_yTL)
% imshow(handles.noSymbol)
% axis image


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CorAlgoGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CorAlgoGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in polyMake.
function polyMake_Callback(hObject, eventdata, handles)
% hObject    handle to polyMake (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


cla(handles.check_yTL)
cla(handles.check_xTL)
cla(handles.check_yBL)
cla(handles.check_xBL)
cla(handles.check_yBR)
cla(handles.check_xBR)
cla(handles.check_yTR)
cla(handles.check_xTR)

set(handles.check_yTL,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');
set(handles.check_xTL,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');
set(handles.check_yBL,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');
set(handles.check_xBL,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');
set(handles.check_yBR,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');
set(handles.check_xBR,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');
set(handles.check_yTR,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');
set(handles.check_xTR,'XTickLabel',[],'YTickLabel',[],...
    'XTick',[],'YTick',[],'Visible','off');

set(handles.yVal_TL,'String','0')
set(handles.xVal_TL,'String','0')
set(handles.yVal_BL,'String','0')
set(handles.xVal_BL,'String','0')
set(handles.yVal_BR,'String','0')
set(handles.xVal_BR,'String','0')
set(handles.yVal_TR,'String','0')
set(handles.xVal_TR,'String','0')

cla(handles.polyDraw)
axes(handles.polyDraw)
imshow(handles.image2show);


[dim1,dim2,~] = size(handles.image2show);
[~, Xcoords, Ycoords] = roipoly(handles.image2show);
mNtb_mask = poly2mask(Xcoords,Ycoords,dim1,dim2);

hold on
plot(Xcoords, Ycoords,'-r');

[B,~,~,~] = bwboundaries(mNtb_mask);

centroid = regionprops(mNtb_mask,handles.image2show(:,:,1),'Centroid');
hold on
plot(centroid.Centroid(1),centroid.Centroid(2), 'b.','MarkerSize', 20);

xmax = max(Xcoords);

ymax = max(Ycoords);

xCor = Xcoords(2:end);
yCor = Ycoords(2:end);

yCorR = zeros(numel(yCor),1);
for yi = 1:numel(yCor)
    yCorR(yi) = yCor(yi) + rand;
end

yCor = yCorR;

xCorR = zeros(numel(xCor),1);
for xi = 1:numel(xCor)
    xCorR(xi) = xCor(xi) + rand;
end

xCor = xCorR;

% Top Left corner
topLeft = (yCor < ymax*0.5 & xCor < xmax*0.5);
yratio = 0.5;
xratio = 0.5;

if sum(topLeft) ~= 1
    if sum(topLeft) == 0
        
        while sum(topLeft) ~= 1
            yratio = yratio + 0.001;
            yval = ymax*yratio;
            set(handles.yVal_TL,'String',num2str(yval))
            xratio = xratio + 0.001;
            xval = ymax*xratio;
            set(handles.xVal_TL,'String',num2str(xval))
            if sum(yCor < yval) == 0
                topLeft = (yCor == min(yCor) & xCor < xval);
            elseif sum(xCor < xval) == 0
                topLeft = (yCor < yval & xCor == min(xCor));
            else
                topLeft = (yCor < yval & xCor < xval);
            end
            % ADD TO THE REST of the VALUES
            updATA = [round(yCor(topLeft)*100)/100 ,  round(xCor(topLeft)*100)/100];
            
            set(handles.polyTable,'Data',updATA, 'ColumnName', {'YVals','XVals'})
            pause(0.05)
        end
        
    elseif sum(topLeft) > 1 % too liberal
        while sum(topLeft) ~= 1
            
            yratio = yratio - 0.0001;
            if yratio < 0
                yratio = 0.5;
            end
            yval = ymax*yratio;
            set(handles.yVal_TL,'String',num2str(yval))
            xratio = xratio - 0.0001;
            if xratio < 0
                xratio = 0.5;
            end
            xval = xmax*xratio;
            set(handles.xVal_TL,'String',num2str(xval))
            if sum(yCor < yval) == 0
                topLeft = (yCor == min(yCor) & xCor < xval);
            elseif sum(xCor < xval) == 0
                topLeft = (yCor < yval & xCor == min(xCor));
            else
                topLeft = (yCor < yval & xCor < xval);
            end
            
            updATA = [round(yCor(topLeft)*100)/100 ,  round(xCor(topLeft)*100)/100];
            
            set(handles.polyTable,'Data',updATA, 'ColumnName', {'YVals','XVals'})
            pause(0.05)
        end
    end
end

axes(handles.check_yTL)
imshow(handles.yesSymbol)
axis image

axes(handles.check_xTL)
imshow(handles.yesSymbol)
axis image

axes(handles.polyDraw)
hold on
plot(xCor(topLeft),yCor(topLeft), 'y.', 'MarkerSize', 30);

% Bottom Left corner
botLeft = (yCor > ymax*0.5 & xCor < xmax*0.5);
yratio = 0.5;
xratio = 0.5;

if sum(botLeft) ~= 1
    if sum(botLeft) == 0 % too conservative
        while sum(botLeft) ~= 1
            yratio = yratio - 0.001;
            yval = ymax*yratio;
            xratio = xratio + 0.001;
            xval = ymax*xratio;
            set(handles.yVal_BL,'String',num2str(yval))
            set(handles.xVal_BL,'String',num2str(xval))
            if sum(yCor > yval) == 0
                botLeft = (yCor == max(yCor) & xCor < xval);
            elseif sum(xCor < xval) == 0
                botLeft = (yCor > yval & xCor == min(xCor));
            else
                botLeft = (yCor > yval & xCor < xval);
            end
            updATA = [round(yCor(botLeft)*100)/100 ,  round(xCor(botLeft)*100)/100];
            
            set(handles.polyTable,'Data',updATA, 'ColumnName', {'YVals','XVals'})
            pause(0.05)
        end
    elseif sum(botLeft) > 1 % too liberal
        while sum(botLeft) ~= 1
            yratio = yratio + 0.001;
            yval = ymax*yratio;
            xratio = xratio - 0.001;
            xval = ymax*xratio;
            set(handles.yVal_BL,'String',num2str(yval))
            set(handles.xVal_BL,'String',num2str(xval))
            if sum(yCor > yval) == 0
                botLeft = (yCor == max(yCor) & xCor < xval);
            elseif sum(xCor < xval) == 0
                botLeft = (yCor > yval & xCor == min(xCor));
            else
                botLeft = (yCor > yval & xCor < xval);
            end
            updATA = [round(yCor(botLeft)*100)/100 ,  round(xCor(botLeft)*100)/100];
            
            set(handles.polyTable,'Data',updATA, 'ColumnName', {'YVals','XVals'})
            pause(0.05)
        end
    end
end

axes(handles.check_yBL)
imshow(handles.yesSymbol)
axis image

axes(handles.check_xBL)
imshow(handles.yesSymbol)
axis image

axes(handles.polyDraw)
hold on
plot(xCor(botLeft),yCor(botLeft), 'y.', 'MarkerSize', 30);

% Bottom Right corner
botRight = (yCor > ymax*0.5 & xCor > xmax*0.5);
yratio = 0.5;
xratio = 0.5;

brCount = 1;
if sum(botRight) ~= 1
    if sum(botRight) == 0
        while sum(botRight) ~= 1 % too conservative
            yratio = yratio - 0.001;
            xratio = xratio - 0.001;
            yval = ymax*yratio;
            xval = xmax*xratio;
            set(handles.yVal_BR,'String',num2str(yval))
            set(handles.xVal_BR,'String',num2str(xval))
            botRight = (yCor > yval & xCor > xval);
            brCount = brCount + 1;
            if brCount > 1000;
                warndlg('bottom right is an infinite loop')
                axes(handles.check_yBR)
                imshow(handles.noSymbol)
                axis image
                break
            end
            updATA = [round(yCor(botRight)*100)/100 ,  round(xCor(botRight)*100)/100];
            
            set(handles.polyTable,'Data',updATA, 'ColumnName', {'YVals','XVals'})
            pause(0.05)
        end
    elseif sum(botRight) > 1
        while sum(botRight) ~= 1 % too liberal
            yratio = yratio + 0.0001;
            yval = ymax*yratio;
            xval = xmax*xratio;
            
            
            botRight = (yCor > yval & xCor > xval);
            botRight = xCor == max(xCor(botRight)) & yCor > yval;
            
            set(handles.yVal_BR,'String',num2str(yval))
            set(handles.xVal_BR,'String',num2str(xval))
            
            brCount = brCount + 1;
            if brCount > 1000;
                warndlg('bottom right is an infinite loop')
                axes(handles.check_yBR)
                imshow(handles.noSymbol)
                axis image
                break
            end
            updATA = [round(yCor(botRight)*100)/100 ,  round(xCor(botRight)*100)/100];
            
            set(handles.polyTable,'Data',updATA, 'ColumnName', {'YVals','XVals'})
            pause(0.05)
        end
    end
end

axes(handles.check_yBR)
imshow(handles.yesSymbol)
axis image

axes(handles.check_xBR)
imshow(handles.yesSymbol)
axis image

axes(handles.polyDraw)
hold on
plot(xCor(botRight),yCor(botRight), 'y.', 'MarkerSize', 30);

% Top Right corner
topRight = (yCor < ymax*0.5 & xCor > xmax*0.5);
yratio = 0.5;
xratio = 0.5;

trCount = 1;
if sum(topRight) ~= 1
    if sum(topRight) == 0
        
        if sum(yCor < ymax*0.5) == 0
            
            maxSortx = sort(xCor,'descend');
            tri = 1;
            while sum(topRight) ~= 1 % too conservative
                maxNow = maxSortx(1:tri);
                xVals = ismember(xCor, maxNow);
                yratio = yratio + 0.001;
                yval = ymax*yratio;
                set(handles.yVal_TR,'String',num2str(yval))
                topRight = (xVals & yCor < yval);
                tri = 1 + 1;
                
                if trCount > 1000;
                    warndlg('bottom right is an infinite loop')
                    axes(handles.check_yTR)
                    imshow(handles.noSymbol)
                    axis image
                    break
                end
                trCount = trCount + 1;
                
            updATA = [round(yCor(topRight)*100)/100 ,  round(xCor(topRight)*100)/100];
            
            set(handles.polyTable,'Data',updATA, 'ColumnName', {'YVals','XVals'})
            pause(0.05)
                
            end
            
        else
            
            while sum(topRight) ~= 1 % too conservative
                yratio = yratio + 0.001;
                xratio = xratio - 0.001;
                
                yval = ymax*yratio;
                xval = xmax*xratio;
                
                set(handles.yVal_TR,'String',num2str(yval))
                set(handles.xVal_TR,'String',num2str(xval))
                
                topRight = (yCor < yval & xCor > xval);
                
                if trCount > 1000;
                    warndlg('bottom right is an infinite loop')
                    axes(handles.check_yTR)
                    imshow(handles.noSymbol)
                    axis image
                    break
                end
                trCount = trCount + 1;
                
            updATA = [round(yCor(topRight)*100)/100 ,  round(xCor(topRight)*100)/100];
            
            set(handles.polyTable,'Data',updATA, 'ColumnName', {'YVals','XVals'})
            pause(0.05)
                
                
            end
            
        end
    elseif sum(topRight) > 1
        
        while sum(topRight) ~= 1 % too liberal
            yratio = yratio - 0.001;
            xratio = xratio + 0.001;
            
            yval = ymax*yratio;
            xval = xmax*xratio;
            
            set(handles.yVal_TR,'String',num2str(yval))
            set(handles.xVal_TR,'String',num2str(xval))
            
            
            topRight = (yCor < ymax*yratio & xCor > xmax*xratio);
            
            if trCount > 1000;
                warndlg('bottom right is an infinite loop')
                axes(handles.check_yTR)
                imshow(handles.noSymbol)
                axis image
                break
            end
            trCount = trCount + 1;
            updATA = [round(yCor(topRight)*100)/100 ,  round(xCor(topRight)*100)/100];
            
            set(handles.polyTable,'Data',updATA, 'ColumnName', {'YVals','XVals'})
            pause(0.05)
            
            
        end
    end
end

axes(handles.check_yTR)
imshow(handles.yesSymbol)
axis image

axes(handles.check_xTR)
imshow(handles.yesSymbol)
axis image

axes(handles.polyDraw)
hold on
plot(xCor(topRight),yCor(topRight), 'y.', 'MarkerSize', 30);

% Top Coordinate
topMidDist = (xCor(topRight) + xCor(topLeft))/2;
topCoord = find(B{1,1}(:,2) == floor(topMidDist) & B{1,1}(:,1) < ymax*0.6, 1, 'first');

axes(handles.polyDraw)
hold on
plot((B{1,1}(topCoord,2)),(B{1,1}(topCoord,1)), 'g.','MarkerSize', 20);

% Bottom Coordinate
botMidDist = (xCor(botRight) + xCor(botLeft))/2;
botCoord = find(B{1,1}(:,2) == floor(botMidDist) & B{1,1}(:,1) > ymax*0.6,1,'last');

plot((B{1,1}(botCoord,2)),(B{1,1}(botCoord,1)), 'g.','MarkerSize', 20);

% Left Coordinate
leftMidDist = (yCor(topLeft) + yCor(botLeft))/2;
leftCoord = find(B{1,1}(:,1) == floor(leftMidDist) & B{1,1}(:,2) < xmax*0.6,1,'first');

plot((B{1,1}(leftCoord,2)),(B{1,1}(leftCoord,1)), 'g.','MarkerSize', 20);

% Right Coordinate
rightMidDist = (yCor(topRight) + yCor(botRight))/2;
rightCoord = find(B{1,1}(:,1) == floor(rightMidDist) & B{1,1}(:,2) > xmax*0.6, 1, 'last');

plot((B{1,1}(rightCoord,2)),(B{1,1}(rightCoord,1)), 'g.','MarkerSize', 20);


% Set the three anchor vertices (top, center , left)
q1_coord1 = [(B{1,1}(topCoord,2));(B{1,1}(topCoord,1))];
q1_coord2 = [centroid.Centroid(1) ; centroid.Centroid(2)];
q1_coord3 = [(B{1,1}(leftCoord,2));(B{1,1}(leftCoord,1))];

% Create X and Y coordinates out of vertices
quadrant_1_Xcoords = [q1_coord1(1);q1_coord2(1);q1_coord3(1)];
quadrant_1_Ycoords = [q1_coord1(2);q1_coord2(2);q1_coord3(2)];

% If Center is greater than left then look for y's that are above both
% left and center
if quadrant_1_Ycoords(2) > quadrant_1_Ycoords(3)
    q1yIndex = (yCor < quadrant_1_Ycoords(2) & yCor < quadrant_1_Ycoords(3));
else % if left is higher then look above the max
    q1yIndex = yCor < max([quadrant_1_Ycoords(2) ; quadrant_1_Ycoords(3)]);
end

if quadrant_1_Xcoords(1) > quadrant_1_Xcoords(2)
    q1xIndex = (xCor < quadrant_1_Xcoords(1) & xCor < quadrant_1_Xcoords(2));
else % if left is higher then look above the max
    q1xIndex = xCor < min([quadrant_1_Xcoords(1) ; quadrant_1_Xcoords(2)]);
end

q1_test = q1yIndex & q1xIndex;

% Get indices for found coordinates
Quad1_coord_indices = find(q1_test == 1);
Quad1Flip = flipud(Quad1_coord_indices);

% Insert found vertices into quadrant coordinates

quad1Xout = zeros(length(Quad1Flip),1);
quad1Yout = zeros(length(Quad1Flip),1);
for q1i2 = 1:length(Quad1Flip);
    quad1Xout(q1i2) = xCor(Quad1Flip(q1i2));
    quad1Yout(q1i2) = yCor(Quad1Flip(q1i2));
end

quad1Xfinal = [quadrant_1_Xcoords ; quad1Xout ; quadrant_1_Xcoords(1)];
quad1Yfinal = [quadrant_1_Ycoords ; quad1Yout ; quadrant_1_Ycoords(1)];

orderIndex.TL = convhull(quad1Xfinal,quad1Yfinal,'simplify',false);

quadCoords.TL = [quad1Xfinal quad1Yfinal];

samplequad = roipoly(dim2,dim1,quadCoords.TL(orderIndex.TL,2),quadCoords.TL(orderIndex.TL,1));

[Qi, ~] = bwboundaries(samplequad,'noholes');
quadIndices = cell2mat(Qi);

plot(quadIndices(:,1),quadIndices(:,2),'m');


