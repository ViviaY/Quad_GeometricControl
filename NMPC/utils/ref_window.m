function Xref_win = ref_window(Xref_full, t_idx, H)
    % NX = size(Xref_full,1); 
    T = size(Xref_full,2);
    i_end = min(t_idx+H, T);
    Xref_win = Xref_full(:, t_idx:i_end);
    % if size(Xref_win,2) < H+1
    %     last_col = Xref_win(:, end);
    %     Xref_win = [Xref_win, repmat(last_col, 1, H+1-size(Xref_win,2))];
    % end
end