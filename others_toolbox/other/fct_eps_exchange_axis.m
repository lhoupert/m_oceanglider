function [] = fct_eps_exchange_axis(fn1, fn2, ichange)
%==========================================================================
%function [] = fct_eps_exchange_axis(fn1, fn2, ichange)
% 26.05.2012, C. Brandt, San Diego (UCSD, CER)
% Last Change: 28.05.2012: preallocation of variables
%--------------------------------------------------------------------------
% FCT_EPS_EXCHANGE_AXIS is a workaround solution for the problem Matlab has
% using multiple colormaps in one figure.
% This function is a quick and dirty solution. It works without using the
% packages 'freezeColor' or 'cbfreeze'. However, I do not know whether it 
% works with every GhostScript version.
% I call it a dirty solution, because the actual problem Matlab exhibits is
% not solved. Instead the resulting figure is obtained from two exported 
% eps-files of the same figure with different colormaps.
%
% The way it works: From two eps files created from the same figure
% using different colormaps this function reads from both the desired axes 
% and puts it altogether to one eps-file 'out.eps'.
%
% If the figure should have more than two different colormaps the 'out.eps'
% can be again used with a third eps-file of the figure using another
% colormap, and so on.
%--------------------------------------------------------------------------
% INPUT
%   fn1: string of 1st eps-file
%   fn2: string of 2nd eps-file
%   ichange: axis that should be exchanged (from file#2 to file#1)
%            the axis numbers are given by its order when put to the figure
%            (colorbars count as axis too!)
% OUTPUT
%   file: out.eps: is equal to eps-file #1, but diagrams 'ichange' are 
%   replaced with diagrams 'ichange' from eps-file #2.
%--------------------------------------------------------------------------
% EXAMPLE:
% Use two different colormaps in one figure (also exchange the colorbars).
% Figure with two axis each having a colorbar.
% First axis should have colormap hot, second axis colormap jet.
%
% figure; clf
% subplot(2,1,1)
% contourf(peaks(64))
% colorbar
% 
% subplot(2,1,2)
% contourf(peaks(64))
% colorbar
% 
% colormap jet; print('-depsc2', '1.eps')
% colormap hot; print('-depsc2', '2.eps')
%
% % Exchange first axis and its colorbar:
% fct_eps_exchange_axis('1.eps', '2.eps', [1, 2])
% disp('>> out.eps << is created containing the merged figures.')
%==========================================================================


%==========================================================================
% Determine Length of both eps-files (For pre-allocation of variables)
%--------------------------------------------------------------------------
% EPS-file #1
fid = fopen(fn1,'rt');
nLines = 0;
while (fgets(fid) ~= -1),
  nLines = nLines+1;
end
fclose(fid);
LF1 = nLines;

% EPS-file #2
fid = fopen(fn2,'rt');
nLines = 0;
while (fgets(fid) ~= -1),
  nLines = nLines+1;
end
fclose(fid);
LF2 = nLines;
%==========================================================================


%==========================================================================
% Read the two EPS files and extract color information
%--------------------------------------------------------------------------
% Read File #1
fid1 = fopen(fn1, 'rt');
lnr      = 0;
ctr_bdef = 0;
ctr_6w   = 0;
l1_bdef  = nan(LF1,1);
l1_6w    = nan(LF1,1);
t1(LF1).line  = [];
while feof(fid1) == 0
  lnr = lnr + 1;
  tline = fgetl(fid1);
  t1(lnr).line  = tline;

  % Look for "/c" and at end "bdef": Indicator in eps for color definitions
  matches=~isempty(findstr(tline,'bdef')) & ~isempty(findstr(tline, '/c'));
  if matches
    ctr_bdef = ctr_bdef+1;
    l1_bdef(ctr_bdef) = lnr;
  end

  % Look for "6 w": '6 w' is an eps indicator of starting a new axis
  % Number of '6 w' is number of all axis +1
  matches = findstr(tline, '6 w');
  if ~isempty(matches)
    ctr_6w = ctr_6w+1;
    l1_6w(ctr_6w) = lnr;
  end

end
fclose(fid1);
% Store Variables and delete NaN numbers
le1 = lnr;
ind = ~isnan(l1_bdef); l1_bdef = l1_bdef(ind);
ind = ~isnan(l1_6w);   l1_6w   = l1_6w(ind);


% Read File #2
fid2 = fopen(fn2, 'rt');
lnr = 0;
ctr_bdef = 0;
ctr_6w   = 0;
l2_bdef  = nan(LF2,1);
l2_6w    = nan(LF2,1);
t2(LF2).line  = [];
while feof(fid2) == 0
  lnr = lnr + 1;
  tline = fgetl(fid2);
  t2(lnr).line  = tline;

  % Look for "/c" and at end "bdef": Indicator in eps for color definitions
  matches=~isempty(findstr(tline,'bdef')) & ~isempty(findstr(tline, '/c'));
  if matches
    ctr_bdef = ctr_bdef+1;
    l2_bdef(ctr_bdef) = lnr;
  end

  % Look for "6 w": '6 w' is an eps indicator of starting a new axis
  % Number of '6 w' is number of all axis +1
  matches = findstr(tline, '6 w');
  if ~isempty(matches)
    ctr_6w = ctr_6w+1;
    l2_6w(ctr_6w) = lnr;
  end
end
fclose(fid2);
% Store Variables and delete NaN numbers
le2 = lnr;
ind = ~isnan(l2_bdef); l2_bdef = l2_bdef(ind);
ind = ~isnan(l2_6w);   l2_6w   = l2_6w(ind);

%==========================================================================
% Extract the color lines lines l1_col and color numbers c1vec (and 2nd)
%--------------------------------------------------------------------------
% Go through File #1 and search for the colors
l1_colnodef.i = nan(LF1,1);
l1_colnodef.n = nan(LF1,1);
l1_col        = nan(LF1,1);
c1vec         = nan(LF1,1);
ctr_col = 0;
ctr_colnodef = 0;
for i=1:le1
    % Look for "6 w": number of axis +1
    matches = findstr(t1(i).line, 'c');
    
    if ~isempty(matches)
      
      if strcmp(t1(i).line(1),'c') && length(t1(i).line)>1
        num = str2double(t1(i).line(2:end));

        if ~isempty(num) && ~isnan(num)
          % Found a color: c<num>
          ctr_col = ctr_col+1;
          l1_col(ctr_col) = i;
          c1vec( ctr_col) = num;
          
          % Look whether line before shows a definition, if not store line
          ind = l1_bdef==i-1;
          if sum(ind)~=1
            ctr_colnodef = ctr_colnodef+1;
            l1_colnodef.i(ctr_colnodef) = i;
            l1_colnodef.n(ctr_colnodef) = num;
          end

        end

      end

    end    
end
% Delete NaN numbers
  % ind = ~isnan(l1_col); l1_col = l1_col(ind);
  % ind = ~isnan(c1vec);  c1vec  = c1vec(ind);
ind = ~isnan(l1_colnodef.i); l1_colnodef.i = l1_colnodef.i(ind);
ind = ~isnan(l1_colnodef.n); l1_colnodef.n = l1_colnodef.n(ind);


% Go through File #2 and search for the colors
c2vec  = nan(LF2,1);
l2_col = nan(LF2,1);
l2_colnodef.i = nan(LF2,1);
l2_colnodef.n = nan(LF2,1);
ctr_col = 0;
ctr_colnodef = 0;
for i=1:le2
    % Look for c<num>
    matches = findstr(t2(i).line, 'c');
    if ~isempty(matches)
      if strcmp(t2(i).line(1),'c') && length(t2(i).line)>1
        num = str2double(t2(i).line(2:end));

        if ~isempty(num) && ~isnan(num)
          % Found a color: c<num>
          ctr_col = ctr_col + 1;
          l2_col(ctr_col) = i;
          c2vec(ctr_col) = num;

          % Look whether line before shows a definition, if not store line
          ind = l2_bdef==i-1;
          if sum(ind)~=1
            ctr_colnodef = ctr_colnodef+1;
            l2_colnodef.i(ctr_colnodef) = i;
            l2_colnodef.n(ctr_colnodef) = num;
          end

        end

      end
    end    
  
end
% Delete NaN numbers
  % ind = ~isnan(l2_col); l2_col = l2_col(ind);
  % ind = ~isnan(c2vec);  c2vec  = c2vec(ind);
ind = ~isnan(l2_colnodef.i); l2_colnodef.i = l2_colnodef.i(ind);
ind = ~isnan(l2_colnodef.n); l2_colnodef.n = l2_colnodef.n(ind);
%==========================================================================



%==========================================================================
% Go through File #1 and replace repetitions of colors with definition
%--------------------------------------------------------------------------
for i=1:length(l1_colnodef.i)
  L = numel(t1);

% Insert the definition in between: previous line and c<num>
  t1( l1_colnodef.i(i)+1:L+1) = t1( l1_colnodef.i(i):L );
  t1( l1_colnodef.i(i) ).line = t1( l1_bdef( l1_colnodef.n(i)+1 ) ).line;

  % +1 for l1_6w included
  ind = l1_6w >= l1_colnodef.i(i);
  l1_6w(ind) = l1_6w(ind) +1;

  % +1 for l1_bdef included
  ind = l1_bdef > l1_colnodef.i(i);
  l1_bdef(ind) = l1_bdef(ind) +1;

  % +1 for l1_colnodef.i included
  ind = l1_colnodef.i >= l1_colnodef.i(i);
  l1_colnodef.i(ind) = l1_colnodef.i(ind) +1;
end

% Go through File #2 and replace repetitions of colors with definition
%--------------------------------------------------------------------------
% For loop: Go through all not defined lines: l2_colnodef
for i=1:length(l2_colnodef.i)
  L = numel(t2);
  
% Insert the definition in between: previous line and c<num>
  t2( l2_colnodef.i(i)+1:L+1) = t2(l2_colnodef.i(i):L);
  t2( l2_colnodef.i(i) ).line = t2( l2_bdef( l2_colnodef.n(i)+1 ) ).line;
  
  % +1 for l2_6w included
  ind = l2_6w >= l2_colnodef.i(i);
  l2_6w(ind) = l2_6w(ind) +1;

  % +1 for l2_bdef included
  ind = l2_bdef > l2_colnodef.i(i);
  l2_bdef(ind) = l2_bdef(ind) +1;

  % +1 for l2_colnodef.i included
  ind = l2_colnodef.i >= l2_colnodef.i(i);
  l2_colnodef.i(ind) = l2_colnodef.i(ind) +1;
  
end
%==========================================================================



%==========================================================================
% Exchange axis/diagrams between the "axis indicator lines" '6 w': ichange
%--------------------------------------------------------------------------
for ic=ichange
% Store t1 in helptext
  ht = t1; clear t1;
  L_ht = numel(ht);
% Position of Picture #1
  a0 = l1_6w(ic  ) +1;
  a1 = l1_6w(ic+1) -1;
  daend = L_ht-(a1+1);
% Position of Picture #2
  b0 = l2_6w(ic  ) +1;
  b1 = l2_6w(ic+1) -1;
  db = b1-b0;

% Insert Begin
  t1(1:a0-1) = ht(1:a0-1);
% Insert Diagramm from File #2
  t1(a0:a0+db) = t2(b0:b1);
% Insert Rest from File #1 stored in helptext
  t1( a0+db+1 : a0+db+1 + daend ) = ht( a1+1 : L_ht );

% Correct '6 w' position due to exchanged figure
  ind = l1_6w > l1_6w(ic);
  l1_6w(ind) = l1_6w(ind) + (b1-b0)-(a1-a0);
end
%==========================================================================



%--------------------------------------write destination file
fid3 = fopen('out.eps', 'wt');
  for i= 1:length(t1)
      fprintf(fid3, '%s\n', t1(i).line );
  end
fclose(fid3);

end