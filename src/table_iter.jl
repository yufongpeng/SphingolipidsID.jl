replace_table!(ex) =  
    @match ex begin
        Expr(:call, fn, Expr(:(->), arg, body), Expr(:(::), tbl, :Table)) => Expr(:call, fn, Expr(:(->), 
                                                                            :index_of_table, 
                                                                            replace_arg!(body, arg, tbl)
                                                                            ), 
                                                                        Expr(:call, :eachindex, tbl))
        
        Expr(:do, Expr(:call, fn, Expr(:(::), tbl, :Table)), Expr(:(->), Expr(:tuple, arg), body)) => Expr(
            :do, Expr(:call, fn, Expr(:call, :eachindex, tbl)), Expr(:(->), Expr(:tuple, :index_of_table), 
            replace_arg!(body, arg, tbl))
        )
        Expr(head, args...) => Expr(head, map(replace_table!, args)...)
        e => e
    end

replace_arg!(body, arg, tbl) = 
    @match body begin
        Expr(:(.), &arg, QuoteNode(at)) => Expr(:ref, Expr(:(.), tbl, QuoteNode(at)), :index_of_table)
        Expr(head, args...) => Expr(head, map(ex -> replace_arg!(ex, arg, tbl), args)...)
        e => e
    end

macro table_iter(fn)
    fn = @match fn begin
        Expr(:(=), Expr(:call, args...), bodies) => Expr(:(=), Expr(:call, args...), bodies)
        Expr(:function, Expr(:call, args...), bodies) => Expr(:function, Expr(:call, args...), bodies)
        Expr(:(=), f, Expr(:(->), args, bodies)) => Expr(:(=), f, Expr(:(->), args, bodies))
    end
    return quote $fn end
end