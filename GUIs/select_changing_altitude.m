function varargout = select_changing_altitude(varargin)
% SELECT_CHANGING_ALTITUDE: A dialog box that allows the user to select
% data ranges from a Merge file where the
%      SELECT_CHANGING_ALTITUDE, by itself, creates a new SELECT_CHANGING_ALTITUDE or raises the existing
%      singleton*.
%
%      H = SELECT_CHANGING_ALTITUDE returns the handle to a new SELECT_CHANGING_ALTITUDE or the handle to
%      the existing singleton*.
%
%      SELECT_CHANGING_ALTITUDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECT_CHANGING_ALTITUDE.M with the given input arguments.
%
%      SELECT_CHANGING_ALTITUDE('Property','Value',...) creates a new SELECT_CHANGING_ALTITUDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before select_changing_altitude_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to select_changing_altitude_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help select_changing_altitude

% Last Modified by GUIDE v2.5 14-Jul-2014 09:51:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @select_changing_altitude_OpeningFcn, ...
                   'gui_OutputFcn',  @select_changing_altitude_OutputFcn, ...
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


% --- Executes just before select_changing_altitude is made visible.
function select_changing_altitude_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to select_changing_altitude (see VARARGIN)

% Choose default command line output for select_changing_altitude
handles.output = hObject;

% If a Merge file is passed to the wrapping function, open the gui with
% that file loaded.
if numel(varargin)>0;
handles = load_merge(handles,varargin{1});
end

% Turn off the buttons not needed until a range is created
% set(handles.set_range_start,'Enable','off');
% set(handles.set_range_end, 'Enable','off');
% set(handles.start_val_text, 'Enable','off');
% set(handles.end_val_text, 'Enable', 'off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes select_changing_altitude wait for user response (see UIRESUME)
uiwait(handles.alt_range_gui);


% --- Outputs from this function are returned to the command line.
function varargout = select_changing_altitude_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Output the handle for the GUI as a second argument
if isfield(handles,'file') && isfield(handles,'ranges')
    varargout{1} = handles.ranges;
    varargout{2} = handles.file;
    varargout{3} = handles.output;
else
    varargout{1} = handles.output;
end
delete(hObject);



% --- Executes on button press in load_merge.
function load_merge_Callback(hObject, eventdata, handles)
% hObject    handle to load_merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[merge_file, merge_path] = uigetfile('.mat','Select a Merge file to open');
if merge_file ~= 0
    handles = load_merge(handles, fullfile(merge_path,merge_file));
    guidata(hObject,handles)
else
    set(handles.feedback_text,'String','Cancelled loading merge file');
end

% --- Executes on button press in set_range_start.
function set_range_start_Callback(hObject, eventdata, handles)
% hObject    handle to set_range_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.range_index)
    ri = handles.range_index;
    set(handles.feedback_text,'String','Select a new starting point for this range');
    [x,y] = ginput(1);
    handles.ranges(ri,1) = x;
    set(handles.start_val_text,'String',sprintf('%.0f',handles.ranges(ri,1)));
    handles = redraw_ranges(handles,x);
    set(handles.feedback_text,'String','');
end
guidata(hObject,handles);

% --- Executes on button press in set_range_end.
function set_range_end_Callback(hObject, eventdata, handles)
% hObject    handle to set_range_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.range_index)
    ri = handles.range_index;
    set(handles.feedback_text,'String','Select a new ending point for this range');
    [x,y] = ginput(1);
    handles.ranges(ri,2) = x;
    set(handles.end_val_text,'String',sprintf('%.0f',handles.ranges(ri,2)));
    handles = redraw_ranges(handles,x);
    set(handles.feedback_text,'String','');
end
guidata(hObject,handles);

% --- Executes on button press in delete_range.
function delete_range_Callback(hObject, eventdata, handles)
% hObject    handle to delete_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.range_index)
    ri = handles.range_index;
    set(handles.feedback_text,'String',sprintf('Deleted range %.0f to %.0f',handles.ranges(ri,1), handles.ranges(ri,2)));
    handles.ranges(ri,:) = [];
    handles = redraw_ranges(handles);
    set(handles.start_val_text,'String','');
    set(handles.end_val_text,'String','');
    handles.range_index = 0;
end
guidata(hObject,handles);


function start_val_text_Callback(hObject, eventdata, handles)
% hObject    handle to start_val_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_val_text as text
%        str2double(get(hObject,'String')) returns contents of start_val_text as a double
if ~isempty(handles.range_index);   
    ri = handles.range_index;
    start_val = str2double(get(hObject,'String'));
    set(handles.feedback_text,'String',sprintf('Changed start value from %.0f to %.0f',handles.ranges(ri,1),start_val));
    handles.ranges(ri,1) = start_val;
    handles = redraw_ranges(handles, start_val);
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function start_val_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_val_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','0');



function end_val_text_Callback(hObject, eventdata, handles)
% hObject    handle to end_val_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_val_text as text
%        str2double(get(hObject,'String')) returns contents of end_val_text as a double
if ~isempty(handles.range_index);   
    ri = handles.range_index;
    end_val = str2double(get(hObject,'String'));
    set(handles.feedback_text,'String',sprintf('Changed end value from %.0f to %.0f',handles.ranges(ri,2),end_val));
    handles.ranges(ri,2) = end_val;
    handles = redraw_ranges(handles, end_val);
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function end_val_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_val_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',0);


% --- Executes on button press in clear_ranges.
function clear_ranges_Callback(hObject, eventdata, handles)
% hObject    handle to clear_ranges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = clearRanges(handles);
guidata(hObject,handles);

% --- Executes on button press in new_range.
function new_range_Callback(hObject, eventdata, handles)
% hObject    handle to new_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.feedback_text,'String','Click on the axes to choose the START range.');
[x,y] = ginput(1);
new_range(1) = x;
set(handles.start_val_text,'String',sprintf('%.0f',x));

set(handles.feedback_text,'String','Click on the axes to choose the END range.');
[x,y] = ginput(1);
new_range(2) = x;
set(handles.end_val_text,'String',sprintf('%.0f',x));
set(handles.feedback_text,'String','');
% 
% set(handles.set_range_start,'Enable','on');
% set(handles.set_range_end, 'Enable','on');
% set(handles.start_val_text, 'Enable','on');
% set(handles.end_val_text, 'Enable', 'on');

handles.ranges = [handles.ranges; new_range];
handles = redraw_ranges(handles,mean(new_range));
handles.range_index = size(handles.ranges,1);
guidata(hObject,handles)


% --- Executes on button press in return_ranges.
function return_ranges_Callback(hObject, eventdata, handles)
% hObject    handle to return_ranges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

alt_range_gui_CloseRequestFcn(handles.alt_range_gui,eventdata,handles);

% --- Executes when user attempts to close alt_range_gui.
function alt_range_gui_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to alt_range_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(handles.alt_range_gui,'waitstatus'),'waiting')
uiresume(handles.alt_range_gui);
else
delete(hObject);
end


% --- Executes on mouse press over axes background.
function alt_plot_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to alt_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pt = get(handles.alt_plot,'CurrentPoint');

% Refresh the drawing of ranges
handles = redraw_ranges(handles, pt);

% Find the range that was just clicked on
if ~isempty(handles.ranges) ri = find(pt(1,1) >= handles.ranges(:,1) & pt(1,1) <= handles.ranges(:,2)); 
else ri = []; 
end

if ~isempty(ri);
    handles.range_index = ri;
    % Set the range text fields to display the currently selected range
    set(handles.start_val_text,'String',sprintf('%.0f',handles.ranges(ri,1)));
    set(handles.end_val_text,'String',sprintf('%.0f',handles.ranges(ri,2)));

%     set(handles.set_range_start,'Enable','on');
%     set(handles.set_range_end, 'Enable','on');
%     set(handles.start_val_text, 'Enable','on');
%     set(handles.end_val_text, 'Enable', 'on');
else
    handles.range_index = ri;
    set(handles.start_val_text,'String','');
    set(handles.end_val_text,'String','');

%     set(handles.set_range_start,'Enable','off');
%     set(handles.set_range_end, 'Enable','off');
%     set(handles.start_val_text, 'Enable','off');
%     set(handles.end_val_text, 'Enable', 'off');
end

guidata(hObject,handles);

% --- Executes on selection change in AltField_popup.
function AltField_popup_Callback(hObject, eventdata, handles)
% hObject    handle to AltField_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AltField_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AltField_popup


% --- Executes during object creation, after setting all properties.
function AltField_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AltField_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','ALTP');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%          Additional functions         %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = redraw_ranges(S, varargin)
if nargin > 1; 
    pt = varargin{1}; 
    pt_x = pt(1,1); 
end
% Remove all existing fills
if ~isempty(S.fills);
    delete(S.fills)
    S.fills = [];
end

ranges = S.ranges;
plot_ylim = S.ylims;
for a=1:size(ranges,1)
    if nargin > 1 && pt_x >= ranges(a,1) && pt_x <= ranges(a,2)
        fillcol = 'b';
    else
        fillcol = 'r';
    end
    fill_x = [ranges(a,1), ranges(a,1), ranges(a,2), ranges(a,2)];
    fill_y = [plot_ylim(1), plot_ylim(2), plot_ylim(2), plot_ylim(1)];
    fh = fill(fill_x, fill_y, fillcol,'FaceAlpha',0.4,'HitTest','off');
    S.fills = [S.fills, fh];
end



function S = clearRanges(S)
if ~isempty(S.fills);
    delete(S.fills)
    S.fills = [];
end
handles.range_index = 0;

% set(handles.set_range_start,'Enable','off');
%     set(handles.set_range_end, 'Enable','off');
%     set(handles.start_val_text, 'Enable','off');
%     set(handles.end_val_text, 'Enable', 'off');

S.ranges = [];

function S = load_merge(S, filename)
load(filename,'Merge');
[fieldname,ok] = select_alt_field(Merge);
if ok
    S.file = filename;
    utc = Merge.Data.UTC.Values;
    alt = eval(sprintf('Merge.Data.%s.Values',fieldname));
    hold off
    alt_line = line(utc,alt,'Parent',S.alt_plot,'color',[0 0.5 0],'HitTest','off');
    hold on
    S.lines.alt_line = alt_line;
    S.ylims = get(S.alt_plot,'YLim');
    if ~isfield(S,'fills'); S.fills = []; end
    S = clearRanges(S);
    set(S.feedback_text,'String',sprintf('%s loaded',filename));
    S.range_index = 0;
else
    set(S.feedback_text,'String','Canceled loading merge file.');
end

function [fieldname, ok] = select_alt_field(Merge)
all_fields = fieldnames(Merge.Data);
if isfield(Merge.Data,'ALTP')
    xx = find(strcmp('ALTP',all_fields));
else
    xx = 1;
end
[fieldindex,ok] = listdlg('ListString',all_fields,'SelectionMode','single','PromptString','Select the altitude field','InitialValue',xx);
fieldname = all_fields{fieldindex};