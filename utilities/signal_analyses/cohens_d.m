function d_value = cohens_d(x1,x2)
    m_1 = nanmean(x1);
    m_2 = nanmean(x2);
    n_1 = sum(1-isnan(x1));
    n_2 = sum(1-isnan(x2));
    s_1 = nanstd(x1);
    s_2 = nanstd(x2);
    
    d_value = abs(m_1 - m_2) / sqrt( ( (n_1-1)*s_1^2 + (n_2-1)*s_2^2 ) / (n_1+n_2-2) );
    
end