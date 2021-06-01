%function [stdevs, bound, F, avgNeighb] = calcCRBGivenLocsAny(method, x, y, blinds, total, channelParam, measMade, d_thr)
%
%OUTPUTS:
% stdevs:    Bounds for the standard deviation of localization error (sqrt of
%               bound on x-variance + bound on y-variance).
% bound:     The matrix bound on covariance
% F:         The fisher information matrix
% avgNeighb: Average number of neighbors per sensor.  Useful if measMade
%               radius is given for neighbor inclusion matrix.
%
%INPUTS:
% method:   Method must be 'T' (TOA), 'A' (AOA), 'R' (RSS), or 'Q' (QRSS).
% x and y:  Actual coordinate vectors.  The first 'blinds' elements are
%              blindfolded, and the remaining are 'reference' or 'anchor'
%              nodes.  x and y must be row vectors.
% blinds:   Total # of blindfolded devices.  The first 'blinds' devices must
%              be correspond to the blindfolded devices, the rest are references.
% total:    Total # of devices.
% channelParam: Equal to sigmadB/n_p for the case of RSS/QRSS measurements, the dB
%              standard deviation of fading divided by the path loss exponent.
%              Or, for the case of TOA measurements , equal to sigma_d, 
%              ie., the std. dev. of distance measurement error (meters).  For the case of 
%              AOA channelParam = sigma_a, the std. dev. of angle measurement
%              error (radians).  The value of channelParam can be scalar if the channel 
%              parameter is the same for all links, or matrix (total by
%              total) if each link has a different channel parameter
% measMade: Measurement matrix.  measMade(i,j)=1 if i and j make a
%              measurement, or =0 if not.  Default is all ones.  If a
%              scalar is sent for measMade, it is considered to be a
%              radius: if ||z_i-z_j|| < radius, then i and j makes measurements,
%              otherwise they don't.
% d_thr:    Distance thresholds used in QRSS.  Not needed for other
%              methods.  The i'th value indicates the distance at which the mean
%              received power seperates the i'th RSS quantization level and
%              the (i+1)'th quantization level.  Note, length(d_thr) = #levels - 1.
%              Also, the d_thr values must be a row vector, and in order (either 
%              decreasing or increasing).
%
% AUTHOR:  Neal Patwari <npatwari@ece.utah.edu>, April 2004
%          AOA part added by Josh Ash <ashj@ece.osu.edu> 7/28/04
%
% EXAMPLE: >> n=20; x=rand(1,n); y=rand(1,n);
%          >> [s, B, F, a] = calcCRBGivenLocsAny('T', x, y, 15, n, 0.2, 0.75);
%          >> calcCRBGivenLocsAny('R', x, y, 15, n, 2.1);
%          >> calcCRBGivenLocsAny('Q', x, y, 15, n, 2.1, ones(n), 0.8:-0.2:0.2)
%
function [stdevs, bound, F, avgNeighb] = calcCRBGivenLocsAny(method, x, y, blinds, total, channelParam, measMade, d_thr)

method = upper(method);
if (method ~='T') & (method ~='R') & (method ~= 'Q' & (method ~= 'A')),
    error('Method must be T (TOA), A (AOA), R (RSS), or Q (QRSS)');
end
% 1. For incomplete measurements, use symmetric matrix measMade to indicate which pairs
%    made measurements (1 for a measurement, 0 for no measurement).  The default 
%    is that all pairs made measurements. 
%       Or, if measMade is just a scalar distance, consider it to be a range, all
%    devices within this range radius 'make measurements' and those outside of 
%    the range don't.
if ~exist('measMade'),
    measMade = ones(total);
end
if length(measMade) == 1,   % Of course total > 1
    rangeRadiusSqr = measMade^2;
    measMade = zeros(total);
else
    rangeRadiusSqr = -1;  % Key for don't use range radius
end
% 2. For the case of identical channel variation on every link, 
%    accept a scalar value for 'channelParam'.
if length(channelParam) == 1,
    channelParam = ones(total).*channelParam;
end
%    Otherwise, use a channel parameter matrix (only the lower triangle is
%    used)
if (method == 'T') | (method == 'A'),
    sigmaConst = channelParam;  % standard deviation of measured distance or angle error
elseif (method == 'R') | (method == 'Q'),
    sigmaConst = channelParam*(log(10)/10);  % 1/sqrt(b), where 'b' is the constant in [1]
end
if method == 'Q',
    K = length(d_thr)+1;  % The number of quantization levels
    % Make sure the d_thr is in decreasing order
    if d_thr(end) > d_thr(1),
        d_thr = fliplr(d_thr);
    end
end

% 3. Calculate the non-diagonal elements of each of the four sub-blocks of F.
%    Each matrix is a superset of the elements needed in F11, F12, and F22,
%    the additional terms are needed for the next step.  
for k = 2:total,
   el = 1:(k-1);
   deltax = x(k) - x(el);
   deltay = y(k) - y(el);

   dSqr   = deltax.^2 + deltay.^2;
   if method == 'T',
       denom(k,el) = dSqr;
   else % for either AOA, RSS or QRSS
       denom(k,el) = dSqr.^2;
   end
   if rangeRadiusSqr > 0,
       measMade(k,el) = (dSqr <= rangeRadiusSqr);
   end
   frontTerm(k,el) = measMade(k,el) ./ (sigmaConst(k,el).^2);
   if method == 'Q',   % There is an additional multiplicative term for QRSS
       d = sqrt(dSqr);
       for j=el;
           sqrtbLogTerms = log(d(j)./d_thr)./sigmaConst(k,j);
           PhiTerms      = [0, N(sqrtbLogTerms), 1];
           deltaPhi      = PhiTerms(2:K+1)-PhiTerms(1:K);
           goodInd       = find(abs(deltaPhi) > eps);
           ExpTerms      = [0, exp(-0.5.*sqrtbLogTerms.^2), 0];
           deltaExpSqr   = (ExpTerms(2:K+1)-ExpTerms(1:K)).^2;
           frontTerm(k,j)= frontTerm(k,j) * sum(deltaExpSqr(goodInd) ./ deltaPhi(goodInd)) / (2*pi);
       end
   end
   if method == 'A'
       term12(k,el) = frontTerm(k,el) .* (-deltay).*(deltax) ./ denom(k,el);
       term12(el,k) = term12(k,el)';
       term11(k,el) = frontTerm(k,el) .* deltay.*deltay ./ denom(k,el);
       term11(el,k) = term11(k,el)';
       term22(k,el) = frontTerm(k,el) .* deltax.*deltax ./ denom(k,el);
       term22(el,k) = term22(k,el)';
   else
       term12(k,el) = frontTerm(k,el) .* deltax.*deltay ./ denom(k,el);
       term12(el,k) = term12(k,el)';
       term11(k,el) = frontTerm(k,el) .* deltax.*deltax ./ denom(k,el);
       term11(el,k) = term11(k,el)';
       term22(k,el) = frontTerm(k,el) .* deltay.*deltay ./ denom(k,el);
       term22(el,k) = term22(k,el)';
   end

end
% 4. Calculate the diagonal elements, which are sums of the elements on each column.
v_xx = sum(term11(:, 1:blinds));
v_xy = sum(term12(:, 1:blinds));
v_yy = sum(term22(:, 1:blinds));

% 5. Combine the two, using only the terms for the blind devices.
F11 = -term11(1:blinds, 1:blinds) + diag(v_xx);
F12 = -term12(1:blinds, 1:blinds) + diag(v_xy);
F22 = -term22(1:blinds, 1:blinds) + diag(v_yy);
F = [F11, F12; F12', F22];
bound = inv(F);

% 6. The location estimate stdev bound: sqrt( var(x) + var(y) ).
stdevs = sqrt(diag(bound(1:blinds,1:blinds)) + diag(bound(1+blinds:2*blinds,1+blinds:2*blinds)))';

% 7. Avg number of neighbors:  note that only the lower triangle is
%    calculated when using a radius, and self-measurement isn't allowed,
%    but the matrix should be symmetric.
avgNeighb = sum(sum(tril(measMade,-1)))*2/total;