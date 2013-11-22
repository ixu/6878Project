var numClusters = 12,
    clusters = {},
	w = 925,
	h = 550,
	margin = 30,
	startDate = 1, 
	endDate = 400,
	minExpression = 0,
	maxExpression = 1000,
	y = d3.scale.linear().domain([maxExpression, minExpression]).range([0 + margin, h - margin]),
	x = d3.scale.linear().domain([startDate, endDate]).range([0 + margin -5, w]),
	timepoints = [1.0,4.0,21.0,116.0,185.0,186.0,255.0,289.0,290.0,292.0,294.0,297.0,301.0,307.0,311.0,322.0,329.0,369.0,380.0,400];

var color = d3.scale.category20();
var sheet = document.styleSheets[0];
for (var i = 0; i < numClusters; i++) {
    clusters['C' + i] = 'Cluster_' + i;
    d3.select("#filters")
        .append('a')
        .attr('id', 'C' + i)
        .html('Cluster ' + i);
    sheet.addRule("#filters a.C" + i, "background: " + color(i) + "; color: #fff;", 1);
    sheet.addRule(".C" + i+".highlight", "stroke: " + color(i) + ";", 1);
}

var vis = d3.select("#vis")
    .append("svg:svg")
    .attr("width", w)
    .attr("height", h)
    .append("svg:g");
    // .attr("transform", "translate(0, 600)");
			
var line = d3.svg.line()
    .x(function(d,i) { return x(d.x); })
    .y(function(d) { return y(d.y); });
					
var genes_clusters = {};
for (var k = 0; k < numClusters; k++) {
    getCluster(k);
}

vis.append("svg:line")
    .attr("x1", x(startDate))
    .attr("y1", y(minExpression))
    .attr("x2", x(endDate))
    .attr("y2", y(minExpression))
    .attr("class", "axis")

vis.append("svg:line")
    .attr("x1", x(startDate))
    .attr("y1", y(minExpression))
    .attr("x2", x(startDate))
    .attr("y2", y(maxExpression))
    .attr("class", "axis")
			
vis.selectAll(".xLabel")
    .data(x.ticks(5))
    .enter().append("svg:text")
    .attr("class", "xLabel")
    .text(String)
    .attr("x", function(d) { return x(d) })
    .attr("y", h-10)
    .attr("text-anchor", "middle")

vis.selectAll(".yLabel")
    .data(y.ticks(4))
    .enter().append("svg:text")
    .attr("class", "yLabel")
    .text(String)
	.attr("x", 0)
	.attr("y", function(d) { return y(d) })
	.attr("text-anchor", "right")
	.attr("dy", 3)
			
vis.selectAll(".xTicks")
    .data(x.ticks(5))
    .enter().append("svg:line")
    .attr("class", "xTicks")
    .attr("x1", function(d) { return x(d); })
    .attr("y1", y(minExpression))
    .attr("x2", function(d) { return x(d); })
    .attr("y2", y(minExpression)+7)
	
vis.selectAll(".yTicks")
    .data(y.ticks(4))
    .enter().append("svg:line")
    .attr("class", "yTicks")
    .attr("y1", function(d) { return y(d); })
    .attr("x1", x(0.5))
    .attr("y2", function(d) { return y(d); })
    .attr("x2", x(startDate))

function onclick(d, i) {
    var currClass = d3.select(this).attr("class");
    if (d3.select(this).classed('selected')) {
        d3.select(this).attr("class", currClass.substring(0, currClass.length-9));
    } else {
        d3.select(this).classed('selected', true);
    }
}

function onmouseover(d, i) {
    var currClass = d3.select(this).attr("class");
    d3.select(this)
        .attr("class", currClass + " current");
    
    var geneName = $(this).attr("name");
    var geneDescription = 'this is a gene';
    //var geneDescription = descriptions[geneName];
    var blurb = '<h2>' + geneName + '</h2>';
    blurb += '<p>' + geneDescription + '</p>'
    
    $("#default-blurb").hide();
    $("#blurb-content").html(blurb);
}
function onmouseout(d, i) {
    var currClass = d3.select(this).attr("class");
    var prevClass = currClass.substring(0, currClass.length-8);
    d3.select(this)
        .attr("class", prevClass);
    $("#default-blurb").show();
    $("#blurb-content").html('');
}

function showCluster(clusterCode) {
    var genes = d3.selectAll("path."+clusterCode);
    if (genes.classed('highlight')) {
        genes.attr("class", clusterCode);
    } else {
        genes.classed('highlight', true);
    }
}

function getCluster(num) {
    d3.text('ProjectCode_v1/Output/Cluster_' + num + '.csv', 'text/csv', function(text) {
    var genes = d3.csv.parseRows(text);
    var average = new Array(20);
    for (var a = 0; a < 20; a++) {
        average[a] = 0;
    }
    for (i=1; i < genes.length; i++) {
        var values = genes[i].slice(1, genes[i.length-1]);
        var currData = [];
        genes_clusters[genes[i][0]] = 'C' + num;

        for (j=0; j < values.length; j++) {
            if (values[j] != '') {
                currData.push({ x: timepoints[j], y: values[j] });  
                average[j] += values[j]/genes.length;
            }
        }
        vis.append("svg:path")
            .data([currData])
            .attr("name", genes[i][0])
            .attr("class", genes_clusters[genes[i][0]])
            .attr("d", line)
            .on("mouseover", onmouseover)
            .on("mouseout", onmouseout);
    }
    for (l = 0; l< average.length; l++) {
        average[l] = {x: timepoints[l], y: average[l]};
    }
    vis.append("svg:path")
    .data([average])
    .attr("name", 'Cluster ' + num + ' Average')
    .attr("class", 'average C' + num)
    .attr("d", line)
    .on("mouseover", onmouseover)
    .on("mouseout", onmouseout);

}); 
}
