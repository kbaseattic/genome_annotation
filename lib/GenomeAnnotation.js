

function GenomeAnnotation(url) {

    var _url = url;


    this.genomeTO_to_reconstructionTO = function(genomeTO)
    {
	var resp = json_call_ajax_sync("GenomeAnnotation.genomeTO_to_reconstructionTO", [genomeTO]);
//	var resp = json_call_sync("GenomeAnnotation.genomeTO_to_reconstructionTO", [genomeTO]);
        return resp[0];
    }

    this.genomeTO_to_reconstructionTO_async = function(genomeTO, _callback, _error_callback)
    {
	json_call_ajax_async("GenomeAnnotation.genomeTO_to_reconstructionTO", [genomeTO], 1, _callback, _error_callback)
    }

    this.genomeTO_to_feature_data = function(genomeTO)
    {
	var resp = json_call_ajax_sync("GenomeAnnotation.genomeTO_to_feature_data", [genomeTO]);
//	var resp = json_call_sync("GenomeAnnotation.genomeTO_to_feature_data", [genomeTO]);
        return resp[0];
    }

    this.genomeTO_to_feature_data_async = function(genomeTO, _callback, _error_callback)
    {
	json_call_ajax_async("GenomeAnnotation.genomeTO_to_feature_data", [genomeTO], 1, _callback, _error_callback)
    }

    this.reconstructionTO_to_roles = function(reconstructionTO)
    {
	var resp = json_call_ajax_sync("GenomeAnnotation.reconstructionTO_to_roles", [reconstructionTO]);
//	var resp = json_call_sync("GenomeAnnotation.reconstructionTO_to_roles", [reconstructionTO]);
        return resp[0];
    }

    this.reconstructionTO_to_roles_async = function(reconstructionTO, _callback, _error_callback)
    {
	json_call_ajax_async("GenomeAnnotation.reconstructionTO_to_roles", [reconstructionTO], 1, _callback, _error_callback)
    }

    this.reconstructionTO_to_subsystems = function(reconstructionTO)
    {
	var resp = json_call_ajax_sync("GenomeAnnotation.reconstructionTO_to_subsystems", [reconstructionTO]);
//	var resp = json_call_sync("GenomeAnnotation.reconstructionTO_to_subsystems", [reconstructionTO]);
        return resp[0];
    }

    this.reconstructionTO_to_subsystems_async = function(reconstructionTO, _callback, _error_callback)
    {
	json_call_ajax_async("GenomeAnnotation.reconstructionTO_to_subsystems", [reconstructionTO], 1, _callback, _error_callback)
    }

    this.annotate_genome = function(genomeTO)
    {
	var resp = json_call_ajax_sync("GenomeAnnotation.annotate_genome", [genomeTO]);
//	var resp = json_call_sync("GenomeAnnotation.annotate_genome", [genomeTO]);
        return resp[0];
    }

    this.annotate_genome_async = function(genomeTO, _callback, _error_callback)
    {
	json_call_ajax_async("GenomeAnnotation.annotate_genome", [genomeTO], 1, _callback, _error_callback)
    }

    this.call_RNAs = function(genomeTO)
    {
	var resp = json_call_ajax_sync("GenomeAnnotation.call_RNAs", [genomeTO]);
//	var resp = json_call_sync("GenomeAnnotation.call_RNAs", [genomeTO]);
        return resp[0];
    }

    this.call_RNAs_async = function(genomeTO, _callback, _error_callback)
    {
	json_call_ajax_async("GenomeAnnotation.call_RNAs", [genomeTO], 1, _callback, _error_callback)
    }

    this.call_CDSs = function(genomeTO)
    {
	var resp = json_call_ajax_sync("GenomeAnnotation.call_CDSs", [genomeTO]);
//	var resp = json_call_sync("GenomeAnnotation.call_CDSs", [genomeTO]);
        return resp[0];
    }

    this.call_CDSs_async = function(genomeTO, _callback, _error_callback)
    {
	json_call_ajax_async("GenomeAnnotation.call_CDSs", [genomeTO], 1, _callback, _error_callback)
    }

    this.find_close_neighbors = function(genomeTO)
    {
	var resp = json_call_ajax_sync("GenomeAnnotation.find_close_neighbors", [genomeTO]);
//	var resp = json_call_sync("GenomeAnnotation.find_close_neighbors", [genomeTO]);
        return resp[0];
    }

    this.find_close_neighbors_async = function(genomeTO, _callback, _error_callback)
    {
	json_call_ajax_async("GenomeAnnotation.find_close_neighbors", [genomeTO], 1, _callback, _error_callback)
    }

    this.assign_functions_to_CDSs = function(genomeTO)
    {
	var resp = json_call_ajax_sync("GenomeAnnotation.assign_functions_to_CDSs", [genomeTO]);
//	var resp = json_call_sync("GenomeAnnotation.assign_functions_to_CDSs", [genomeTO]);
        return resp[0];
    }

    this.assign_functions_to_CDSs_async = function(genomeTO, _callback, _error_callback)
    {
	json_call_ajax_async("GenomeAnnotation.assign_functions_to_CDSs", [genomeTO], 1, _callback, _error_callback)
    }

    this.annotate_proteins = function(genomeTO)
    {
	var resp = json_call_ajax_sync("GenomeAnnotation.annotate_proteins", [genomeTO]);
//	var resp = json_call_sync("GenomeAnnotation.annotate_proteins", [genomeTO]);
        return resp[0];
    }

    this.annotate_proteins_async = function(genomeTO, _callback, _error_callback)
    {
	json_call_ajax_async("GenomeAnnotation.annotate_proteins", [genomeTO], 1, _callback, _error_callback)
    }

    function _json_call_prepare(url, method, params, async_flag)
    {
	var rpc = { 'params' : params,
		    'method' : method,
		    'version': "1.1",
	};
	
	var body = JSON.stringify(rpc);
	
	var http = new XMLHttpRequest();
	
	http.open("POST", url, async_flag);
	
	//Send the proper header information along with the request
	http.setRequestHeader("Content-type", "application/json");
	//http.setRequestHeader("Content-length", body.length);
	//http.setRequestHeader("Connection", "close");
	return [http, body];
    }

    /*
     * JSON call using jQuery method.
     */

    function json_call_ajax_sync(method, params)
    {
        var rpc = { 'params' : params,
                    'method' : method,
                    'version': "1.1",
        };
        
        var body = JSON.stringify(rpc);
        var resp_txt;
	var code;
        
        var x = jQuery.ajax({       "async": false,
                                    dataType: "text",
                                    url: _url,
                                    success: function (data, status, xhr) { resp_txt = data; code = xhr.status },
				    error: function(xhr, textStatus, errorThrown) { resp_txt = xhr.responseText, code = xhr.status },
                                    data: body,
                                    processData: false,
                                    type: 'POST',
				    });

        var result;

        if (resp_txt)
        {
	    var resp = JSON.parse(resp_txt);
	    
	    if (code >= 500)
	    {
		throw resp.error;
	    }
	    else
	    {
		return resp.result;
	    }
        }
	else
	{
	    return null;
	}
    }

    function json_call_ajax_async(method, params, num_rets, callback, error_callback)
    {
        var rpc = { 'params' : params,
                    'method' : method,
                    'version': "1.1",
        };
        
        var body = JSON.stringify(rpc);
        var resp_txt;
	var code;
        
        var x = jQuery.ajax({       "async": true,
                                    dataType: "text",
                                    url: _url,
                                    success: function (data, status, xhr)
				{
				    resp = JSON.parse(data);
				    var result = resp["result"];
				    if (num_rets == 1)
				    {
					callback(result[0]);
				    }
				    else
				    {
					callback(result);
				    }
				    
				},
				    error: function(xhr, textStatus, errorThrown)
				{
				    if (xhr.responseText)
				    {
					resp = JSON.parse(xhr.responseText);
					if (error_callback)
					{
					    error_callback(resp.error);
					}
					else
					{
					    throw resp.error;
					}
				    }
				},
                                    data: body,
                                    processData: false,
                                    type: 'POST',
				    });

    }

    function json_call_async(method, params, num_rets, callback)
    {
	var tup = _json_call_prepare(_url, method, params, true);
	var http = tup[0];
	var body = tup[1];
	
	http.onreadystatechange = function() {
	    if (http.readyState == 4 && http.status == 200) {
		var resp_txt = http.responseText;
		var resp = JSON.parse(resp_txt);
		var result = resp["result"];
		if (num_rets == 1)
		{
		    callback(result[0]);
		}
		else
		{
		    callback(result);
		}
	    }
	}
	
	http.send(body);
	
    }
    
    function json_call_sync(method, params)
    {
	var tup = _json_call_prepare(url, method, params, false);
	var http = tup[0];
	var body = tup[1];
	
	http.send(body);
	
	var resp_txt = http.responseText;
	
	var resp = JSON.parse(resp_txt);
	var result = resp["result"];
	    
	return result;
    }
}

