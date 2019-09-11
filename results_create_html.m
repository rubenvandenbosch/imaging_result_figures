function results_create_html(fileName,title,imgFiles,borders)
% create_html(imgFiles)
% 
% DESCRIPTION
% Create HTML file with collection of all figures in input imgFiles
% 
% INPUT
% - fileName : char; File name of output HTML file
% - title    : structure with title for html file and headings printed at 
%              the top.
% - imgFile  : cellstr; list of image files to include in the html file
% - borders  : struct with fields "big" and "small". Contains search string 
%              for image file paths to base separation of images with big 
%              and small borders on. Set both to false for no borders.
% 
% -------------------------------------------------------------------------
% Ruben van den Bosch
% Donders Institute, Radboud University Nijmegen
% July 2019
% 

% Create HTML file and write content
% -------------------------------------------------------------------------
% Create file
fid = fopen(fileName,'w+');

% Opening stuff
fprintf(fid,'<!DOCTYPE html>\n<html>\n<head>\n');

% Title
HTMLtitle = title.title;
fprintf(fid,sprintf('<title>%s</title>\n',HTMLtitle));

% Headings: title
head1 = title.head1;
fprintf(fid,sprintf('<h1>%s</h1>\n',head1));
head2 = title.head2;
fprintf(fid,sprintf('<h2>%s</h2>\n',head2));

% Close head and open style
fprintf(fid,'</head>\n<style>\n');

% Scale images using css
scale = '100%%';
fprintf(fid,sprintf('img.resize {\n\tmax-width:%s;\n\tmax-heigth:%s;\n}\n',scale,scale));

% CSS settings for big and small borders between images
fprintf(fid,'hr.big {\n\tmargin-top: 2em;\n\tmargin-bottom: 2em;\n\tborder: 3px dashed darkred;\n}\n');
fprintf(fid,'hr.small {\n\tmargin-top: 2em;\n\tmargin-bottom: 2em;\n\tborder: 2px dotted darkred;\n}\n');

% Close style and open body of contents
fprintf(fid,'</style>\n<body>\n');

% Link to image files and print borders
% -------------------------------------------------------------------------
% Initial matches for big and small borders search strings
if borders.big
    big     = regexp(imgFiles{1},filesep,'split');
    big     = big{contains(big(1:end-1),borders.big)};
end
if borders.small
    small   = regexp(imgFiles{1},filesep,'split');
    small   = small{contains(small(1:end-1),borders.small)};
end

for i = 1:numel(imgFiles)
    [~,fname,~] = fileparts(imgFiles{i});

    % If under windows, change path separator to unix style
    if strcmp(filesep,'\')
        imgFiles{i} = strrep(imgFiles{i},'\','/');
    end

    % Print borders
    if borders.big
        temp = regexp(imgFiles{i},'/','split');
        temp = temp{contains(temp(1:end-1),borders.big)};
        
        % Add big border if new big border separator instance found
        if ~strcmp(temp,big)
            fprintf(fid,'<hr class="big">\n');
            big = temp;
            bigprinted = true;
        else
            bigprinted = false;
        end
    end
    if borders.small
        temp = regexp(imgFiles{i},'/','split');
        temp = temp{contains(temp(1:end-1),borders.small)};
        
        % Add small border if new small border separator instance found
        if ~strcmp(temp,small) && ~bigprinted
            fprintf(fid,'<hr class="small">\n');
            small = temp;
        elseif bigprinted
            small = temp;
        end
    end
    
    % Print link to image
    fprintf(fid,sprintf('<img class="resize" src="file:/%s" alt="%s">\n',imgFiles{i},fname));    
end

% Closing statements
fprintf(fid,'</body>\n</html>\n');

% Close file
fclose(fid);
end