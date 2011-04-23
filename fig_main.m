function varargout = fig_main(varargin)
% FIG_MAIN M-file for fig_main.fig
%      FIG_MAIN, by itself, creates a new FIG_MAIN or raises the existing
%      singleton*.
%
%      H = FIG_MAIN returns the handle to a new FIG_MAIN or the handle to
%      the existing singleton*.
%
%      FIG_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIG_MAIN.M with the given input arguments.
%
%      FIG_MAIN('Property','Value',...) creates a new FIG_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fig_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fig_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fig_main

% Last Modified by GUIDE v2.5 08-Apr-2011 16:41:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fig_main_OpeningFcn, ...
                   'gui_OutputFcn',  @fig_main_OutputFcn, ...
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


% --- Executes just before fig_main is made visible.
function fig_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fig_main (see VARARGIN)

% Choose default command line output for fig_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fig_main wait for user response (see UIRESUME)
% uiwait(handles.fig_main);
globals;
Hd1 = 0; Hd2 = 0;
fig_main_pictureData = imread('images/Sxema.png');
Img=image(fig_main_pictureData, 'Parent', handles.axes_Pic1);
set(handles.axes_Pic1, 'XTick', []);
set(handles.axes_Pic1, 'YTick', []);
clc

% --- Outputs from this function are returned to the command line.
function varargout = fig_main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function fig_main_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fig_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pb_Calculate.
function pb_Calculate_Callback(hObject, eventdata, handles)
globals;

Fd = 44.2e6;
Td = 1/Fd;
L = fix(Tmod / Td);
t = (1:L)*Td;

%dfi  = 0e6;
%f0 = 14e6 + dfi;
f0 = str2double(get(handles.ed_f0, 'String')) * 1e6;
phi = str2double(get(handles.ed_phi, 'String')) / 180 * pi;
Gdk = ((1:L)>(L/16)) - ((1:L)>(2*L/16)) + ((1:L)>(3*L/16)) - ((1:L)>(4*L/16)) + ((1:L)>(5*L/16)) - ((1:L)>(6*L/16)) + ((1:L)>(7*L/16))- 0.5;
S = Gdk .* cos(2*pi*f0*t + phi);

fp = str2double(get(handles.ed_fp, 'String')) * 1e6;
H = str2double(get(handles.ed_H, 'String')) * 1e6;


load Hd1.mat
%fvtool(handles.axes12, Hd1,'Color','White');
y_Filt = filter(Hd1, S);

plot(handles.axes1, 1e6*t, S, 1e6*t, y_Filt)
xlim(handles.axes1, [min(1e6*t) max(1e6*t)])
grid(handles.axes1, 'on');

plot(handles.axes2,(0:(Fd/L):(Fd*(1-1/L)))/1e6,abs(fft(y_Filt)))
xlim(handles.axes2,[0 (Fd*(1-1/L))/1e6])
grid(handles.axes2, 'on');

q = y_Filt .* sin(2*pi*f0*t);
i = y_Filt .* cos(2*pi*f0*t);
ic = i; % ����� ���������� �������, �� ����������
plot(handles.axes3, i, q)
maxx = max(abs(i));
if max(abs(q))>maxx
    maxx = max(abs(q));
end
xlim(handles.axes3,[-maxx*1.2 +maxx*1.2])
ylim(handles.axes3,[-maxx*1.2 +maxx*1.2])
grid(handles.axes3, 'on');

load Hd2.mat
i_Filt = filter(Hd2, i);
q_Filt = filter(Hd2, q);

plot(handles.axes4, (0:(Fd/L):(Fd*(1-1/L)))/1e6, abs(fft(i_Filt)), (0:(Fd/L):(Fd*(1-1/L)))/1e6, abs(fft(q_Filt)))
xlim(handles.axes4, [0 (Fd*(1-1/L))/1e6])
grid(handles.axes4, 'on');

plot(handles.axes5, i_Filt, q_Filt)
maxx = max(abs(i_Filt));
if max(abs(q_Filt))>maxx
    maxx = max(abs(q_Filt));
end
xlim(handles.axes5, [-maxx*1.2 +maxx*1.2])
ylim(handles.axes5, [-maxx*1.2 +maxx*1.2])
grid(handles.axes5, 'on');

plot(handles.axes6, 1e6*t, rad2deg(atan2(q_Filt,i_Filt))/360*187, 1e6*t, -phi/(2*pi)*187 - (Gdk + 0.5 - 1) * pi/(2*pi)*187)
grid(handles.axes6, 'on');
xlim(handles.axes6, [min(1e6*t) max(1e6*t)])

plot(handles.axes16, 1e6*t, i_Filt, 1e6*t, q_Filt)
grid(handles.axes16, 'on');
xlim(handles.axes16, [min(1e6*t) max(1e6*t)])

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
globals;
hF=figure(1);
plot(1e6*t, S, 1e6*t, y_Filt)
xlim([min(1e6*t) max(1e6*t)])
xlabel('t, {\mu}s')
ylabel('S, y_{Filt}')
grid on

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
globals;
hF=figure(2);
plot((0:(Fd/L):(Fd*(1-1/L)))/1e6,abs(fft(y_Filt)))
xlabel('f, MHz')
ylabel('|fft(y_{Filt})|')
xlim([0 (Fd*(1-1/L))/1e6])
grid on

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
globals;
hF=figure(4);
plot((0:(Fd/L):(Fd*(1-1/L)))/1e6, abs(fft(i_Filt)), (0:(Fd/L):(Fd*(1-1/L)))/1e6, abs(fft(q_Filt)))
xlabel('f, MHz')
ylabel('Okno_{LPH}, fft i, fft q')
xlim([0 (Fd*(1-1/L))/1e6])
grid on

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
globals;
hF=figure(6);
plot(1e6*t, rad2deg(atan2(q_Filt,i_Filt))/360*187, 1e6*t, -phi/(2*pi)*187 - (Gdk + 0.5 - 1) * pi/(2*pi)*187)
xlabel('t, \mu')
ylabel('atan Filt, mm')
grid on

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
globals;
hF=figure(5);
plot(i_Filt, q_Filt)
xlabel('i_{Filt}')
ylabel('q_{Filt}')
maxx = max(abs(i_Filt));
if max(abs(q_Filt))>maxx
    maxx = max(abs(q_Filt));
end
xlim([-maxx*1.2 +maxx*1.2])
ylim([-maxx*1.2 +maxx*1.2])
grid on

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
globals;
hF=figure(3);
plot(ic, q)
xlabel('i')
ylabel('q')
maxx = max(abs(ic));
if max(abs(q))>maxx
    maxx = max(abs(q));
end
xlim([-maxx*1.2 +maxx*1.2])
ylim([-maxx*1.2 +maxx*1.2])
grid on



function ed_phi_Callback(hObject, eventdata, handles)
% hObject    handle to ed_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_phi as text
%        str2double(get(hObject,'String')) returns contents of ed_phi as a double


% --- Executes during object creation, after setting all properties.
function ed_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_f0_Callback(hObject, eventdata, handles)
% hObject    handle to ed_f0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_f0 as text
%        str2double(get(hObject,'String')) returns contents of ed_f0 as a double


% --- Executes during object creation, after setting all properties.
function ed_f0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_f0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_fp_Callback(hObject, eventdata, handles)
% hObject    handle to ed_fp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_fp as text
%        str2double(get(hObject,'String')) returns contents of ed_fp as a double


% --- Executes during object creation, after setting all properties.
function ed_fp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_fp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_H_Callback(hObject, eventdata, handles)
% hObject    handle to ed_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_H as text
%        str2double(get(hObject,'String')) returns contents of ed_H as a double


% --- Executes during object creation, after setting all properties.
function ed_H_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
globals;
load Hd1.mat


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
globals;
if Hd1 ~= 0
    fvtool(Hd1,'Color','White');
else
    pushbutton8_Callback(0, 0, 0);
    fvtool(Hd1,'Color','White');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
globals;
load Hd2.mat


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
globals;
if Hd2 ~= 0
    fvtool(Hd2,'Color','White');
else
    pushbutton10_Callback(0, 0, 0);
    fvtool(Hd2,'Color','White');
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
globals;
figure(7)
plot(1e6*t, i_Filt, 1e6*t, q_Filt)
xlabel('t, {\mu}s')
ylabel('i_{Filt}, q_{Filt}')
grid('on');
xlim([min(1e6*t) max(1e6*t)])
