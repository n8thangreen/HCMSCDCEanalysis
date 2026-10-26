    library(kableExtra)

    format_table <- function(x, caption = NULL, scale_down = TRUE) {
      if (knitr::is_latex_output()) {
        k <- knitr::kable(x, format = "latex", booktabs = TRUE, caption = caption)
        if (scale_down) {
          k <- kableExtra::kable_styling(k, latex_options = "scale_down")
        }
        return(k)
      } else {
        return(knitr::kable(x, caption = caption))
      }
    }

    df <- data.frame(A = 1:5, B = 6:10, C = 11:15, D = 16:20, E = 21:25, F = 26:30, G = 31:35, H = 36:40)
    format_table(df, caption = "Test Table")

<table>
<caption>Test Table</caption>
<thead>
<tr class="header">
<th style="text-align: right;">A</th>
<th style="text-align: right;">B</th>
<th style="text-align: right;">C</th>
<th style="text-align: right;">D</th>
<th style="text-align: right;">E</th>
<th style="text-align: right;">F</th>
<th style="text-align: right;">G</th>
<th style="text-align: right;">H</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">1</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">26</td>
<td style="text-align: right;">31</td>
<td style="text-align: right;">36</td>
</tr>
<tr class="even">
<td style="text-align: right;">2</td>
<td style="text-align: right;">7</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">22</td>
<td style="text-align: right;">27</td>
<td style="text-align: right;">32</td>
<td style="text-align: right;">37</td>
</tr>
<tr class="odd">
<td style="text-align: right;">3</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">18</td>
<td style="text-align: right;">23</td>
<td style="text-align: right;">28</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">38</td>
</tr>
<tr class="even">
<td style="text-align: right;">4</td>
<td style="text-align: right;">9</td>
<td style="text-align: right;">14</td>
<td style="text-align: right;">19</td>
<td style="text-align: right;">24</td>
<td style="text-align: right;">29</td>
<td style="text-align: right;">34</td>
<td style="text-align: right;">39</td>
</tr>
<tr class="odd">
<td style="text-align: right;">5</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">20</td>
<td style="text-align: right;">25</td>
<td style="text-align: right;">30</td>
<td style="text-align: right;">35</td>
<td style="text-align: right;">40</td>
</tr>
</tbody>
</table>

Test Table
