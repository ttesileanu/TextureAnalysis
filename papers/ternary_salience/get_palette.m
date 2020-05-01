function [c, d] = get_palette()
% GET_PALETTE Get default colors used for the paper.
%   c = GET_PALETTE() returns an N x 3 array of RGB colors used in the
%   paper. The number of colors `N` is currently 7.
%
%   The colors were obtained using the tool at coolors.co.
%
%   c, d = GET_PALETTE() also returns a map that can be used to obtain
%   colors by name. The names currently in the dictionary are 'red',
%   'blue', 'orange', 'green', 'gray', 'dark blue', and 'light orange'.

d = containers.Map;
d('blue') = hex2color('016FB9');            % spanish blue
d('red') = hex2color('CC4D4D');             % english vermillion
d('orange') = hex2color('FF9505');          % yellow orange (color wheel)
d('green') = hex2color('419D78');           % paolo veronese green
d('gray') = hex2color('353531');            % jet
d('dark blue') = hex2color('2D3047');       % space adet
d('light orange') = hex2color('FFDBB5');    % light orange

c = [
    d('blue') ;
    d('red') ;
    d('orange') ;
    d('green') ;
    d('gray') ;
    d('dark blue') ;
    d('light orange')
];

end