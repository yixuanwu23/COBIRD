function [intersect, t, u, v, xcoor] = intersect_triangle_ray (...
  orig, dir, vert0, vert1, vert2, varargin)
% SPDX-FileCopyrightText: 2018 Jaroslaw Tuszynski
% SPDX-License-Identifier: BSD-2-Clause
%
% Adapted for COBIRD from "Triangle/Ray Intersection" by Jaroslaw Tuszynski.
% Original source:
% https://www.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection
% Full license text: LICENSES/BSD-2-Clause.txt
%
%INTERSECT_TRIANGLE_RAY Ray/triangle intersection.
%    INTERSECT = cobird.geometry.intersect_triangle_ray(ORIG, DIR, VERT0, VERT1, VERT2)
%      calculates ray/triangle intersections using the algorithm proposed
%      by Möller and Trumbore (1997) in vectorized MATLAB code. The ray
%      starts at ORIG and points along DIR. The triangle is defined by
%      vertices VERT0, VERT1, and VERT2. All input
%      arrays are in Nx3 or 1x3 format, where N is number of triangles or
%      rays.
%
%   [INTERSECT, T, U, V, XCOOR] = cobird.geometry.intersect_triangle_ray(...)
%     Returns:
%     * Intersect - boolean array of length N informing which line and
%                 triangle pair intersect
%     * t   - distance from the ray origin to the intersection point in
%             units of |dir|. Provided only for line/triangle pair that
%             intersect unless 'fullReturn' parameter is true.
%     * u,v - barycentric coordinates of the intersection point
%     * xcoor - Cartesian coordinates of the intersection point
%
%   cobird.geometry.intersect_triangle_ray(...,'param','value',...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * planeType - 'one sided' or 'two sided' (default) - how to treat
%        triangles. In 'one sided' version only intersections in single
%        direction are counted and intersections with back facing
%           triangles are ignored
%    * lineType - 'ray' (default), 'line' or 'segment' - how to treat rays:
%        - 'line' means infinite (on both sides) line;
%        - 'ray' means infinite (on one side) ray coming out of origin;
%        - 'segment' means the segment from ORIG to ORIG + DIR
%    * border - controls border handling:
%        - 'normal'(default) border - triangle is exactly as defined.
%           Intersections with border points can be easily lost due to
%           rounding errors.
%        - 'inclusive' border - triangle is marginally larger.
%           Intersections with border points are always captured but can
%           lead to double counting when working with surfaces.
%        - 'exclusive' border - triangle is marginally smaller.
%           Intersections with border points are not captured and can
%           lead to under-counting when working with surfaces.
%    * eps - (default = 1e-5) controls border size
%    * fullReturn - (default = false) controls returned variables t, u, v,
%        and xcoor. By default in order to save time, not all t, u & v are
%        calculated, only t, u & v for intersections can be expected.
%        fullReturn set to true will force the calculation of them all.
%
% ALGORITHM:
%  Function solves
%        |t|
%    M * |u| = (o-v0)
%        |v|
%  for [t; u; v] where M = [-d, v1-v0, v2-v0]. u,v are barycentric coordinates
%  and t - the distance from the ray origin in |d| units
%  ray/triangle intersect if u>=0, v>=0 and u+v<=1
%
% NOTE:
%  The algorithm is able to solve several types of problems:
%  * many faces / single ray  intersection
%  * one  face  / many   rays intersection
%  * one  face  / one    ray  intersection
%  * many faces / many   rays intersection
%  In order to allow that to happen all input arrays are expected in Nx3
%  format, where N is number of vertices or rays. In most cases number of
%  vertices is different than number of rays, so one of the inputs will
%  have to be cloned to have the right size. Use "repmat(A,size(B,1),1)".
%
% Based on:
%  * "Fast, Minimum Storage Ray-Triangle Intersection," Tomas Möller and
%    Ben Trumbore, Journal of Graphics Tools, 2(1):21--28, 1997.
%    http://www.graphics.cornell.edu/pubs/1997/MT97.pdf
%  * http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/raytri/
%  * http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/raytri/raytri.c
%
% Original author: Jaroslaw (Jarek) Tuszynski
% Copyright (c) 2018, Jaroslaw Tuszynski. All rights reserved.
% Licensed under the BSD 2-Clause License.

%% Transpose inputs if needed
if (size(orig ,1)==3 && size(orig ,2)~=3), orig =orig' ; end
if (size(dir  ,1)==3 && size(dir  ,2)~=3), dir  =dir'  ; end
if (size(vert0,1)==3 && size(vert0,2)~=3), vert0=vert0'; end
if (size(vert1,1)==3 && size(vert1,2)~=3), vert1=vert1'; end
if (size(vert2,1)==3 && size(vert2,2)~=3), vert2=vert2'; end

%% In case of single points clone them to the same size as the rest
N = max([size(orig,1), size(dir,1), size(vert0,1), size(vert1,1), size(vert2,1)]);
if (size(orig ,1)==1 && N>1 && size(orig ,2)==3), orig  = repmat(orig , N, 1); end
if (size(dir  ,1)==1 && N>1 && size(dir  ,2)==3), dir   = repmat(dir  , N, 1); end
if (size(vert0,1)==1 && N>1 && size(vert0,2)==3), vert0 = repmat(vert0, N, 1); end
if (size(vert1,1)==1 && N>1 && size(vert1,2)==3), vert1 = repmat(vert1, N, 1); end
if (size(vert2,1)==1 && N>1 && size(vert2,2)==3), vert2 = repmat(vert2, N, 1); end

%% Check if all the sizes match
SameSize = (any(size(orig)==size(vert0)) && ...
  any(size(orig)==size(vert1)) && ...
  any(size(orig)==size(vert2)) && ...
  any(size(orig)==size(dir  )) );
assert(SameSize && size(orig,2)==3, ...
  'All input vectors have to be in Nx3 format.');

%% Read user preferences
eps        = 1e-5;
planeType  = 'two sided';
lineType   = 'ray';
border     = 'normal';
fullReturn = false;
nVarargs   = length(varargin);
assert(mod(nVarargs, 2) == 0, 'Options must be supplied as name-value pairs.');
k = 1;
while (k<=nVarargs)
    assert(ischar(varargin{k}) || (isstring(varargin{k}) && isscalar(varargin{k})), ...
        'Option names must be character vectors or string scalars.')
    switch lower(char(varargin{k}))
      case 'eps'
        eps = abs(varargin{k+1});
        k = k+1;
      case 'planetype'
        planeType = lower(strtrim(varargin{k+1}));
        k = k+1;
      case 'border'
        border = lower(strtrim(varargin{k+1}));
        k = k+1;
      case 'linetype'
        lineType = lower(strtrim(varargin{k+1}));
        k = k+1;
      case 'fullreturn'
        fullReturn = (double(varargin{k+1})~=0);
        k = k+1;
      otherwise
        error('cobird:geometry:intersect_triangle_ray:unknownOption', ...
            'Unknown option: %s.', char(varargin{k}));
    end
    k = k+1;
end

%% Set up border parameter
switch border
  case 'normal'
    zero=0.0;
  case 'inclusive'
    zero=eps;
  case 'exclusive'
    zero=-eps;
  otherwise
    error('Border parameter must be either "normal", "inclusive" or "exclusive"')
end

%% initialize default output
intersect = false(size(orig,1),1); % by default there are no intersections
t = inf+zeros(size(orig,1),1); u=t; v=t;
xcoor = nan+zeros(size(orig));

%% Find faces parallel to the ray
edge1 = vert1-vert0;          % find vectors for two edges sharing vert0
edge2 = vert2-vert0;
tvec  = orig -vert0;          % vector from vert0 to ray origin
pvec  = cross(dir, edge2,2);  % begin calculating determinant - also used to calculate U parameter
det   = sum(edge1.*pvec,2);   % determinant of the matrix M = dot(edge1,pvec)
switch planeType
  case 'two sided'            % treats triangles as two sided
    angleOK = (abs(det)>eps); % if determinant is near zero then ray lies in the plane of the triangle
  case 'one sided'            % treats triangles as one sided
    angleOK = (det>eps);
  otherwise
    error('Triangle parameter must be either "one sided" or "two sided"');
end
if all(~angleOK), return; end % if all parallel than no intersections

%% Different behavior depending on one or two sided triangles
det(~angleOK) = nan;              % change to avoid division by zero
u    = sum(tvec.*pvec,2)./det;    % 1st barycentric coordinate
if fullReturn
  % calculate all variables for all line/triangle pairs
  qvec = cross(tvec, edge1,2);    % prepare to test V parameter
  v    = sum(dir  .*qvec,2)./det; % 2nd barycentric coordinate
  t    = sum(edge2.*qvec,2)./det; % 'position on the line' coordinate
  % test if line/plane intersection is within the triangle
  ok   = (angleOK & u>=-zero & v>=-zero & u+v<=1.0+zero);
else
  % limit some calculations only to line/triangle pairs where it makes
  % a difference. It is tempting to try to push this concept of
  % limiting the number of calculations to only the necessary to "u"
  % and "t" but that produces slower code
  v = nan+zeros(size(u)); t=v;
  ok = (angleOK & u>=-zero & u<=1.0+zero); % mask
  % if all line/plane intersections are outside the triangle than no intersections
  if ~any(ok), intersect = ok; return; end
  qvec = cross(tvec(ok,:), edge1(ok,:),2); % prepare to test V parameter
  v(ok,:) = sum(dir(ok,:).*qvec,2) ./ det(ok,:); % 2nd barycentric coordinate
  if (~strcmpi(lineType,'line')) % 'position on the line' coordinate
    t(ok,:) = sum(edge2(ok,:).*qvec,2)./det(ok,:);
  end
  % test if line/plane intersection is within the triangle
  ok = (ok & v>=-zero & u+v<=1.0+zero);
end

%% Test where along the line the line/plane intersection occurs
switch lineType
  case 'line'      % infinite line
    intersect = ok;
  case 'ray'       % ray is bound on one side
    intersect = (ok & t>=-zero); % intersection on the correct side of the origin
  case 'segment'   % segment is bound on two sides
    intersect = (ok & t>=-zero & t<=1.0+zero); % intersection between origin and destination
  otherwise
    error('lineType parameter must be either "line", "ray" or "segment"');
end

%% calculate intersection coordinates if requested
if (nargout>4)
  ok = intersect | fullReturn;
  xcoor(ok,:) = vert0(ok,:) ...
    + edge1(ok,:).*repmat(u(ok,1),1,3) ...
    + edge2(ok,:).*repmat(v(ok,1),1,3);
end
